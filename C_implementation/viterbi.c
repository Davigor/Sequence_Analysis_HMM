/**
 * Simple implementation of the viterbi algorithm for estimating the states of a Hidden Markov Model given a sequence text file.
 * Program assumes there are 2 states, state 1 and state 2. State transition matrix probabilites and emission lambda for sampling
 * from Poisson distribution can be altered by user.
 * 
 * Optional argument to read in file of known states for comparison with algorithm's output. 
 * Sequence file is assumed to be one entry per line, and state file is assumed to give corresponding state on same line separated
 * whitespace (see .txt files for example).
 * 
 * Usage: ./viterbi n my_sequence_file.txt my_state_file.txt
 * n = number of entries in sequence file (required)
 * my_sequence_file.txt = sequence file (required)
 * my_state_file.txt = state file (optional)
 *
**/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double max (double a, double b);
double poisson (double k, double lambda);
int argmax (double row0, double row1);

int main (int argc, char *argv[]) 
{
    // check for correct number of command line args
    if (argc != 3 && argc != 4)
    {
        printf("Usage: ./viterbi n my_sequence_file.txt my_state_file.txt .  Include at least n and sequence file.\n");
        return 1;
    }
    
    int n = atoi(argv[1]);
    int num;
    
    // open sequence file and store in array
    FILE *seqf = fopen(argv[2], "r");
    
    if (seqf == NULL)
    {
        printf("Invalid sequence file.\n");
        return 1;
    }
    
    int *seq = malloc(n * sizeof(int));
    
    for (int i = 0; i < n; i++)
    {
        fscanf(seqf, "%i", &num);    
        seq[i] = num;
    }
    fclose(seqf);
    
    // if passed as an argument, open the state solution file and store as an array
    if(argv[3])
    {
        FILE *statef = fopen(argv[3], "r");
        
        if (statef == NULL)
        {
            printf("Invalid state file.\n");
            return 1;
        }
        
        int *state = calloc(n, sizeof(int));
        
        for (int i = 0; i < n; i++)
        {
            fscanf(statef, "%*i %i", &num);
            state[i] = num;
            printf("%i", state[i]);
        }
        fclose(statef);
        free(state);
        printf("\n");
    }

    // state transition matrix in log space
    double a[2][2] = {
        { log(0.9551),  log(0.0449) },
        { log(0.0880),  log(0.9120) }
    };
    
    // emission lambdas for sampling from poisson distribution
    double e[2] = {1.8234, 5.7812};
    
    // allocate memory
    int *path = calloc(n, sizeof(double));
    double **vprob = calloc(n, sizeof *vprob);
    double **ptr = calloc(n, sizeof *ptr);
    double **pi = calloc(n, sizeof *pi);
    
    for (int i = 0; i < 2; i++)
    {
        vprob[i] = calloc(n, sizeof *vprob[i]);
        ptr[i] = calloc(n, sizeof *ptr[i]);
        pi[i] = calloc(n, sizeof *pi[i]);
    }

    // initialize vprob array; assumed starting state is state 1
    vprob[0][0] = 1;
    vprob[1][0] = 0;
    double row0;
    double row1;
    
    // viterbi algorithm in log space to avoid underflow. Emission probabilities sampled from poisson distribution
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            row0 = (vprob[0][i - 1] + a[0][j]);
            row1 = (vprob[1][i - 1] + a[1][j]);
            
            vprob[j][i] = log( poisson(seq[i - 1], e[j]) ) + max( row0, row1 );
            ptr[j][i] = argmax( row0, row1 );
            pi[j][i] = max( row0 , row1 );
        }
    }
    free(seq);
    
    // traceback to find most likely path
    path[n - 1] = argmax( pi[0][n - 1], pi[1][n - 1] );
    for (int i = n - 2; i > 0; i--)
    {
        path[i] = ptr[path[i + 1]][i + 1];
    }

    // free remaining memory
    for (int i = 0; i < 2; i++)
    {
        free(vprob[i]);
        free(ptr[i]);
        free(pi[i]);
    }
    free(vprob);
    free(ptr);
    free(pi);

    // incremment each entry of the path to correspond to states 1 and 2 intead of 0 and 1, and print
    for (int i = 0; i < n; i++)
    {
        path[i]++;
        printf("%i", path[i]);
    }
    printf("\n");
    free(path);
    return 0;
}

double max (double a, double b)
{
    if (a > b)
    {
        return a;
    }
    else if (a < b)
    {
    return b;
    }
    // if equal, returns arbitrary argument for specific use in this algorithm
    return b;
}

int argmax (double row0, double row1)
{
    if (row0 > row1)
    {
        return 0;
    }
    else if (row0 < row1)
    {
        return 1;
    }
    return row1;
}

double poisson (double k, double lambda)
{
    int factk;
    if (k == 0)
    {
        factk = 1;
    }
    else
    {
        factk = k;
        for (int i = k - 1; i > 0; i--)
        {
            factk *= i;
        }
    }
    return ( (pow(lambda, k) * exp(-lambda)) / factk );
}
