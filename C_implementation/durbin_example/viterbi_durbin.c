/**
 * Implementation of the viterbi algorithm for estimating the states of a Hidden Markov Model given at least a sequence text file.
 * Program automatically determines n value from sequence file and assumes that state file has same n value.
 * 
 * Program follows example from Durbin et. al.'s "The occasionally dishonest casino, part 1." on p. 54 of Biological Sequence
 * Analyis, with solution and viterbi output given on p. 57. The two states, F and L, correspond to a "Fair" or a "Loaded" die.
 * 
 * free .pdf of Durbin at: http://dnapunctuation.org/~poptsova/course/Durbin-Et-Al-Biological-Sequence-Analysis-CUP-2002-No-OCR.pdf 
 * 
 * Optional argument to read in file of known states for comparison with algorithm's output. 
 * Sequence file and state files are is assumed to be one entry per line (see .txt files for example).
 * 
 * Usage: ./viterbi my_sequence_file.txt [my_state_file.txt]
 *
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <err.h>
#include <sysexits.h>

int sequence_length = 0;
int* sequence;
int* states;

static double max (double a, double b)
{
    if (a > b)
    {
        return a;
    }
    // if equal, returns arbitrary argument for specific use in this algorithm
    else
    {
        return b;
    }
}

static int argmax (double row0, double row1)
{
    if (row0 > row1)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

static int* read_sequencefile(const char* sequence_file, int* out_n)
{
    FILE *seqf = fopen(sequence_file, "r");
    if (!seqf)
    {
        fprintf(stderr, "Invalid sequence file.\n");
        return NULL;
    }
    
    int num;
    int n = 0;
    sequence = malloc(sizeof(int));
    
    while(fscanf(seqf, "%i", &num) == 1)
    {
        sequence[n] = num - 1;
        n++;
        sequence = realloc(sequence, (n + 1) * sizeof(int));
        
        if(!sequence)
        {
            fprintf(stderr, "Not enough memory.");
            return NULL;
        }
    }
    fclose(seqf);
    
    *out_n = n;
    return sequence;
}


static int read_statefile(const char* state_file)
{
    FILE *statef = fopen(state_file, "r");
    
    if (!statef)
    {
        fprintf(stderr, "Invalid state file.\n");
        return EXIT_FAILURE;
    }
    
    char *state = malloc(sequence_length * sizeof(char));
    if(!state)
    {
        fprintf(stderr, "Not enough memory.");
        return EXIT_FAILURE;
    }
    
    char ch;
    printf("State solution:\n");
    for (int i = 0; i < sequence_length; i++)
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
    return 0;
}


static int run_viterbi(int* seq, int seq_length)
{
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
    int *path = malloc(seq_length * sizeof(path));
    double **vprob = malloc(seq_length * sizeof(*vprob));
    double **ptr = malloc(seq_length * sizeof(*ptr));
    double **pi = malloc(seq_length * sizeof(*pi));
    
    if( !path || !vprob || !ptr || !pi )
        {
            fprintf(stderr, "Not enough memory.");
            return 1;
        }
    
    for (int i = 0; i < 2; i++)
    {
        vprob[i] = malloc(seq_length * sizeof(double));
        ptr[i] = malloc(seq_length * sizeof(double));
        pi[i] = malloc(seq_length * sizeof(double));
        
        if( !vprob[i] || !ptr[i] || !pi[i] )
        {
            fprintf(stderr, "Not enough memory.");
            return 1;
        }
    }

    // initialize vprob array; assumed starting state is state F
    vprob[0][0] = 1;
    vprob[1][0] = 0;
    double row0;
    double row1;
    
    // viterbi algorithm in log space to avoid underflow
    for (int i = 1; i < seq_length; i++)
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
    free(sequence);
    
    // traceback to find most likely path
    path[seq_length - 1] = argmax( pi[0][seq_length - 1], pi[1][seq_length - 1] );
    for (int i = seq_length - 2; i > 0; i--)
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
    for (int i = 0; i < seq_length; i++)
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

int main (int argc, char *argv[]) 
{
    // check for correct number of command line args
    if (argc != 2 && argc != 3)
    {
        fprintf(stderr, "Usage: ./viterbi my_sequence_file.txt [my_state_file.txt]. Include at least sequence file.\n");
        return 1;
    }
    
    read_sequencefile(argv[1], &sequence_length);
    
    if(argv[2])
    {
        read_statefile(argv[2]);
    }
    
    run_viterbi(sequence, sequence_length);

    return 0;
}
