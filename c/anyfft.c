/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o anyfft anyfft.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./anyfft
 **************************************************************************************************/

/**************************************************************************************************
 Includes necessary libraries:
 **************************************************************************************************/
#include <stdlib.h>                            // Standard Library;
#include <stdio.h>                             // Input and Output;
#include <math.h>                              // Math Functions;
#include <time.h>                              // Time measurement;

// #include <stdlib.h>                    // Biblioteca padrão da linguagem;
// #include <stdio.h>                     // Entrada e saída;
// #include <sys/types.h>                 // Tipos de dados;
// #include <string.h>                    // Processamento de strings;
// #include <ctype.h>                     // Processamento de caracteres;
// #include <math.h>                      // Funções matemáticas;


/**************************************************************************************************
 Definitions:
 **************************************************************************************************/
#define REPEAT 500                             // Number of executions to compute average time;


/**************************************************************************************************
 Data structure to deal with complex numbers. There is, of course, an external library that
 implements operations for comple numbers, but I'm defining this here to reduce dependencies.
 **************************************************************************************************/
typedef struct {
    float r, i;                                // Real and Imaginary parts of a complex number;
} Complex;



/**************************************************************************************************
 Small set of functions to operate with complex numbers
 **************************************************************************************************/
// Initialization of a complex number:
Complex cmplx(float r, float i)
{
    Complex w;

    w.r = r;                                   // Real part;
    w.i = i;                                   // Imaginary part;
    return w;
}


// Complex addition of numbers a and b:
Complex cadd(Complex a, Complex b)
{
    Complex w;

    w.r = a.r + b.r;                           // Real part of sum;
    w.i = a.i + b.i;                           // Imaginary part of sum;
    return w;
}


// Complex subtraction of number b from a:
Complex csub(Complex a, Complex b)
{
    Complex w;

    w.r = a.r - b.r;                           // Real part of difference;
    w.i = a.i - b.i;                           // Imaginary part of difference;
    return w;
}


// Complex product of numbers a and b:
Complex cmul(Complex a, Complex b)
{
    Complex w;

    w.r = a.r*b.r - a.i*b.i;                   // Real part of product;
    w.i = a.r*b.i + a.i*b.r;                   // Imaginary part of product;
    return w;
}


// Complex exponential of an angle:
Complex cexpn(float angle)
{
    Complex w;

    w.r = cos(angle);                          // Real part of exponential;
    w.i = sin(angle);                          // Imaginary part of exponential;
    return w;
}


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
void complex_show(Complex x[], int n)
{
    for (int i=0; i<n; i++)
        printf("(%7.4f, %7.4f)\n", x[i].r, x[i].i);
}


/**************************************************************************************************
 * Auxiliary function: time_it
 *   This function calls a Fast Fourier Transform function repeatedly a certain number of times,
 *   measure execution time and average it.
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
float time_it(void (*f)(Complex *, Complex *, int), int size, int repeat)
{
    Complex x[1024], X[1024];                  // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++) {                // Initialize the vector;
        x[j].r = j;
        x[j].i = 0;
    }
    t0 = clock();                              // Start counting time;
    for(int j=0; j<repeat; j++)
        (*f)(x, X, size);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeat;
    return t1;
}


/**************************************************************************************************
 * Function: factor
 *   This function finds the smallest prime factor of a given number. If the argument is prime
 *   itself, then it is the return value.
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
    int rn = n/2;                              // Search up to half the number;
    for(int i=2; i<rn; i++)
        if (n%i == 0) return i;                // If remainder is zero, a factor is found;
    return n;
}


/**************************************************************************************************
 * Function: direct_ft
 *   Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that has
 *   O(N^2) complexity.
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
void direct_ft(Complex x[], Complex X[], int N)
{
    Complex W, Wk, Wkn, w;                     // Twiddle factors;

    W = cexpn(-2*M_PI/N);                      // Initializes twiddle factors;
    Wk = cmplx(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = cmplx(0, 0);                    // Accumulates the results;
        Wkn = cmplx(1, 0);                     // Initializes twiddle factors;
        for(int n=0; n<N; n++) {
            w = cmul(x[n], Wkn);
            X[k] = cadd(X[k], w);
            Wkn = cmul(Wkn, Wk);               // Updates twiddle factors;
        }
        Wk = cmul(Wk, W);
    }
}


/**************************************************************************************************
 * Function: recursive_fft
 *   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This has
 *   smallest complexity than the direct FT, though the exact value is difficult to compute.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. It must be a composite number, or else the
 *     computation will be defered to the direct FT, and there will be no efficiency gain.
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **********************************************************************************************/
void recursive_fft(Complex x[], Complex X[], int N)
{
    Complex *xj, *Xj;                          // Subsequences of the transform;
    Complex W, Wj, Wkj;
    int N1, N2;

    N1 = factor(N);                            // Smallest prime factor of the length;
    if (N1==N)                                 // If the length is prime itself,
        direct_ft(x, X, N);                    //   the transform is given by the direct form;
    else {
        N2 = N / N1;                           // Decompose in two factors, N1 being prime;
        xj = malloc(sizeof(Complex)*N2);       // Allocates memory for subsequences
        Xj = malloc(sizeof(Complex)*N2);       //    and their transforms;
        W = cexpn(-2*M_PI/N);                  // Twiddle factor;
        Wj = cmplx(1, 0);
        for(int j=0; j<N1; j++) {              // Computes every subsequence of size N2;
            for(int n=0; n<N2; n++) {
                xj[n] = x[n*N1+j];             // Creates the subsequence;
                Xj[n] = cmplx(0, 0);           // Inicializes the transform;
            }
            recursive_fft(xj, Xj, N2);         // Computes the DFT of the subsequence;
            Wkj = cmplx(1, 0);
            for(int k=0; k<N; k++) {           // Recombine results;
                X[k] = cadd(X[k], cmul(Xj[k%N2], Wkj));
                Wkj = cmul(Wkj, Wj);           // Update twiddle factors;
            }
            Wj = cmul(Wj, W);
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

    // Starts by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | Direta  | Recurs. |\n");
    printf("+---------+---------+---------+---------+\n");

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int i=0; i<8; i++) {

        // Computes the average execution time:
        n = SIZES[i];
        dtime = time_it(direct_ft, n, REPEAT);
        rtime = time_it(recursive_fft, n, REPEAT);

        // Print the results:
        printf("| %7d | %7d | %7.4f | %7.4f |\n", n, n*n, dtime, rtime);
    }
    printf("+---------+---------+---------+---------+\n");

    return 0;
}

