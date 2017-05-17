/**
 * Simple implementation of the viterbi algorithm for estimating the states of a Hidden Markov Model given a sequence text file.
 * Program assumes there are 2 states, state 1 and state 2. State transition matrix probabilites and emission lambda for sampling
 * from Poisson distribution are hardcoded in program. 
 * 
 * Program automatically determines n value from sequence file and assumes that state file has same n value.
 * 
 * Optional argument to read in file of known states for comparison with algorithm's output. 
 * Sequence file is assumed to be one entry per line, and state file is assumed to give corresponding state on same line separated
 * by whitespace (see .txt files for example).
 * 
 * Usage: ./viterbi my_sequence_file.txt my_state_file.txt
 * my_sequence_file.txt = sequence file (required)
 * my_state_file.txt = state file (optional)
 *
**/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

double max (double a, double b);
double poisson (double k, double lambda);
int argmax (double row0, double row1);

int main (int argc, char *argv[]) 
{
    // check for correct number of command line args
    if (argc != 2 && argc != 3)
    {
        printf("Usage: ./viterbi my_sequence_file.txt my_state_file.txt .  Include at least sequence file.\n");
        return 1;
    }
    
    // open sequence file and store in array. Dynamically allocate memory and automatically detect sequence n value
    FILE *seqf = fopen(argv[1], "r");
    if (!seqf)
    {
        printf("Invalid sequence file.\n");
        return 1;
    }
    
    int num;
    int memsize = 100;
    int n = 0;
    int *seq = calloc(memsize, sizeof(int));
    
    while(fscanf(seqf, "%i", &num) == 1)
    {
        seq[n] = num;
        n++;

        if (n == memsize)
        {
            memsize += 100;
            seq = realloc(seq, memsize * sizeof(int));
        }
        
        if(!seq)
        {
            printf("Not enough memory.");
            return 1;
        }
    }
    fclose(seqf);
    
    // if passed as an argument, open the state solution file and print. Assumes n is same as sequence n above
    if(argv[2])
    {
        FILE *statef = fopen(argv[2], "r");
        
        if (!statef)
        {
            printf("Invalid state file.\n");
            return 1;
        }
        
        int *state = calloc(n, sizeof(int));
        if(!state)
        {
            printf("Not enough memory.");
            return 1;
        }
        
        printf("State solution:\n");
        for (int i = 0; i < n; i++)
        {
            fscanf(statef, "%*i" "%i", &num);
            state[i] = num;
            printf("%i", state[i]);
        }
        fclose(statef);
        free(state);
        printf("\n\n");
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
    
    // print most likely path, incrementing by 1 for states 1 and 2
    for (int i = 0; i < n; i++)
    {
        path[i]++;
        printf("%i", path[i]);
    }
    printf("\n");

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
