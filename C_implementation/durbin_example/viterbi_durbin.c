/**
 * Implementation of the viterbi algorithm for estimating the states of a Hidden Markov Model given at least a sequence text file.
 * Program automatically determines n value from sequence file and assumes that state file has same n value.
 * 
 * Program follows example from Durbin et. al.'s "The occasionally dishonest casino, part 1." on p. 54 of Biological Sequence
 * Analyis (2002), with solution and viterbi output given on p. 57. The two states, F and L, correspond to a "Fair" or a "Loaded" die.
 * 
 * free .pdf of Durbin at: http://dnapunctuation.org/~poptsova/course/Durbin-Et-Al-Biological-Sequence-Analysis-CUP-2002-No-OCR.pdf
 * 
 * Optional argument to read in file of known states for comparison with algorithm's output. 
 * Sequence file and state files are is assumed to be one entry per line (see .txt files for example).
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
        seq[n] = num - 1;
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
        
        char *state = calloc(n, sizeof(char));
        if(!state)
        {
            printf("Not enough memory.");
            return 1;
        }
        
        char ch;
        printf("State solution:\n");
        for (int i = 0; i < n; i++)
        {
            fscanf(statef, "%c %*[\r\n]", &ch);
            state[i] = ch;
            
            if (i % 60 == 0 && i != 0)
            {
                printf("\n");
            }
            printf("%c", state[i]);
        }
        fclose(statef);
        free(state);
        printf("\n\n");
    }
    
    // state transition matrix in log space
    double a[2][2] = {
        { log(0.95),  log(0.05) },
        { log(0.1),  log(0.9) }
    };
    
    // emission probabilities, corresponding to p of rolling 1 thru 6 on fair or loaded die
    double e[6][2] = {
        { log( ((double) 1)/6),  log(0.1) },
        { log( ((double) 1)/6),  log(0.1) },
        { log( ((double) 1)/6),  log(0.1) },
        { log( ((double) 1)/6),  log(0.1) },
        { log( ((double) 1)/6),  log(0.1) },
        { log( ((double) 1)/6),  log(0.5) },
    };
    
    // allocate rest of memory and error handle
    int *path = calloc(n, sizeof(double));
    double **vprob = calloc(n, sizeof(double *));
    double **ptr = calloc(n, sizeof(double *));
    double **pi = calloc(n, sizeof(double *));
    
    if( !path || !vprob || !ptr || !pi )
        {
            printf("Not enough memory.");
            return 1;
        }
    
    for (int i = 0; i < 2; i++)
    {
        vprob[i] = calloc(n, sizeof(double));
        ptr[i] = calloc(n, sizeof(double));
        pi[i] = calloc(n, sizeof(double));
        
        if( !vprob[i] || !ptr[i] || !pi[i] )
        {
            printf("Not enough memory.");
            return 1;
        }
    }

    // initialize vprob array; assumed starting state is state F
    vprob[0][0] = 1;
    vprob[1][0] = 0;
    double row0;
    double row1;
    
    // viterbi algorithm in log space to avoid underflow
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            row0 = (vprob[0][i - 1] + a[0][j]);
            row1 = (vprob[1][i - 1] + a[1][j]);
            
            vprob[j][i] = e[seq[i]][j] + max( row0, row1 );
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

    // print viterbi result
    printf("Viterbi output:\n");
    for (int i = 0; i < n; i++)
    {   
        if (i % 60 == 0 && i != 0)
        {
            printf("\n");
        }
        
        if (path[i] == 0)
        {
            printf("F");
        }
        if (path[i] == 1)
        {
            printf("L");
        }
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
