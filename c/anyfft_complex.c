/**************************************************************************************************
 * Fast Fourier Transform -- C with Native Complex Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o anyfft anyfft_complex.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./anyfft
 **************************************************************************************************/

/**************************************************************************************************
 Include necessary libraries:
 **************************************************************************************************/
#include <stdlib.h>                            // Standard Library;
#include <stdio.h>                             // Input and Output;
#include <math.h>                              // Math Functions;
#include <complex.h>                           // Complex numbers;
#include <time.h>                              // Time measurement;


/**************************************************************************************************
 Definitions:
 **************************************************************************************************/
#define REPEAT 500                             // Number of executions to compute average time;


/**************************************************************************************************
 * Auxiliary function: complex_show
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 *   n
 *     Number of elements on the vector.
 **************************************************************************************************/
void complex_show(float complex x[], int n)
{
    for (int i=0; i<n; i++)
        printf("(%7.4f, %7.4f)\n", creal(x[i]), cimag(x[i]));
}


/**************************************************************************************************
 * Auxiliary function: time_it
 *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
 *
 * Parameters:
 *  f
 *    Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 *  size
 *    Number of elements in the vector on which the transform will be applied;
 *  repeat
 *    Number of times the function will be called.
 *
 * Returns:
 *   The average execution time for that function with a vector of the given size.
 **************************************************************************************************/
float time_it(void (*f)(float complex *, float complex *, int), int size, int repeat)
{
    float complex x[1024], X[1024];            // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = j;
    t0 = clock();                              // Start counting time;
    for(int j=0; j<repeat; j++)
        (*f)(x, X, size);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeat;
    return t1;
}


/**************************************************************************************************
 * Function: direct_ft
 *   Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
 *   complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the DFT will be computed. Given the nature of the implementation, there
 *     is no restriction on the size of the vector;
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void direct_ft(float complex x[], float complex X[], int N)
{
    float complex W, Wk, Wkn;                  // Twiddle factors;

    W = cexpf(-2i*M_PI/N);                     // Initialize twiddle factors;
    Wk = 1;
    for(int k=0; k<N; k++) {
        X[k] = 0;                              // Accumulate the results;
        Wkn = 1;                               // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            X[k] = X[k] + Wkn * x[n];
            Wkn = Wkn * Wk;                    // Update twiddle factors;
        }
        Wk = Wk * W;
    }
}


/**************************************************************************************************
 * Function: factor
 *   Smallest prime factor of a given number. If the argument is prime itself, then it is the
 *   return value.
 *
 * Parameters:
 *   n
 *     Number to be inspected.
 *
 * Returns:
 *   The smallest prime factor, or the number itself if it is already a prime.
 **************************************************************************************************/
int factor(int n)
{
    int rn = (int) ceil(sqrt(n));              // Search up to the square root of the number;
    for(int i=2; i<=rn; i++)
        if (n%i == 0) return i;                // If remainder is zero, a factor is found;
    return n;
}


/**************************************************************************************************
 * Function: recursive_fft
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
 *   complexity than the direct FT, though the exact value is difficult to compute.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. Its length must be a composite number, or else
 *     the computation will be defered to the direct FT, and there will be no efficiency gain.
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void recursive_fft(float complex x[], float complex X[], int N)
{
    float complex *xj, *Xj;                    // Subsequences of the transform;
    float complex W, Wj, Wkj;
    int N1, N2;

    N1 = factor(N);                            // Smallest prime factor of the length;
    if (N1==N)                                 // If the length is prime itself,
        direct_ft(x, X, N);                    //   the transform is given by the direct form;
    else {
        N2 = N / N1;                           // Decompose in two factors, N1 being prime;

        xj = malloc(sizeof(float complex)*N2);         // Allocate memory for subsequences
        Xj = malloc(sizeof(float complex)*N2);         //    and their transforms;

        W = cexpf(-2i*M_PI/N);                 // Twiddle factor;
        Wj = 1;
        for(int j=0; j<N1; j++) {              // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++) {
                xj[n] = x[n*N1+j];             // Create the subsequence;
                Xj[n] = 0;                     // Initialize the transform;
            }
            recursive_fft(xj, Xj, N2);         // Compute the DFT of the subsequence;
            Wkj = 1;
            for(int k=0; k<N; k++) {
                X[k] = X[k] + Wkj * Xj[k%N2];  // Recombine results;
                Wkj = Wkj * Wj;                // Update twiddle factors;
            }
            Wj = Wj * W;
        }

        free(xj);                              // Clean-up used memory;
        free(Xj);
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, itime;
    int SIZES[] = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
    int n;

    // Start by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | Direct  | Recurs. |\n");
    printf("+---------+---------+---------+---------+\n");

    // Try it with vectors with the given sizes:
    for(int i=0; i<8; i++) {

        // Compute the average execution time:
        n = SIZES[i];
        dtime = time_it(direct_ft, n, REPEAT);
        rtime = time_it(recursive_fft, n, REPEAT);

        // Print the results:
        printf("| %7d | %7d | %7.4f | %7.4f |\n", n, n*n, dtime, rtime);
    }
    printf("+---------+---------+---------+---------+\n");

    return 0;
}
