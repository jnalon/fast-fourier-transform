/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o fft fft.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fft
 **************************************************************************************************/

/**************************************************************************************************
 Include necessary libraries:
 **************************************************************************************************/
#include <stdlib.h>                            // Standard Library;
#include <stdio.h>                             // Input and Output;
#include <math.h>                              // Math Functions;
#include <time.h>                              // Time measurement;


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
 Small set of functions to operate with complex numbers:
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
float time_it(void (*f)(Complex *, Complex *, int), int size, int repeat)
{
    Complex x[1024], X[1024];                  // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = cmplx(j, 0);
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
 *     is no restriction on the size of the vector, although it will almost always be called with a
 *     power of two size to give a fair comparison;
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void direct_ft(Complex x[], Complex X[], int N)
{
    Complex W, Wk, Wkn, w;                     // Twiddle factors;

    W = cexpn(-2*M_PI/N);                      // Initialize twiddle factors;
    Wk = cmplx(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = cmplx(0, 0);                    // Accumulate the results;
        Wkn = cmplx(1, 0);                     // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            w = cmul(x[n], Wkn);
            X[k] = cadd(X[k], w);
            Wkn = cmul(Wkn, Wk);               // Update twiddle factors;
        }
        Wk = cmul(Wk, W);
    }
}


/**************************************************************************************************
 * Function: recursive_fft
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 *   complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector of
 *     a power of two length, or it will fail. No checks on this are made.
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void recursive_fft(Complex x[], Complex X[], int N)
{
    Complex *xe, *xo, *Xe, *Xo;                // Vectors with intermediate results;
    Complex W, Wk, w;                          // Twiddle factors;
    int N2;

    if(N==1)                                   // A length-1 vector is its own FT;
        X[0] = x[0];
    else {
        N2 = N >> 1;

        xe = (Complex *) malloc(sizeof(Complex)*N2);   // Allocate memory for computation;
        xo = (Complex *) malloc(sizeof(Complex)*N2);
        Xe = (Complex *) malloc(sizeof(Complex)*N2);
        Xo = (Complex *) malloc(sizeof(Complex)*N2);

        for(int k=0; k<N2; k++) {                      // Split even and odd samples;
            xe[k] = x[k<<1];
            xo[k] = x[(k<<1)+1];
        }
        recursive_fft(xe, Xe, N2);                     // Transform of even samples;
        recursive_fft(xo, Xo, N2);                     // Transform of odd samples;

        W = cexpn(-2*M_PI/N);                          // Twiddle factors;
        Wk = cmplx(1, 0);
        for(int k=0; k<N2; k++) {
            w = cmul(Xo[k], Wk);
            X[k] = cadd(Xe[k], w);                     // Recombine results;
            X[k+N2] = csub(Xe[k], w);
            Wk = cmul(Wk, W);                          // Update twiddle factors;
        }

        free(Xo);                                      // Releasing memory of intermediate vectors;
        free(Xe);
        free(xo);
        free(xe);
    }
}


/**************************************************************************************************
 * Function: bit_reverse
 *   Bit-reversed version of an integer number.
 *
 * Parameters:
 *   k
 *     The number to be bit-reversed;
 *   r
 *     The number of bits to take into consideration when reversing.
 *
 * Returns:
 *   The number k, bit-reversed according to integers with r bits.
 **************************************************************************************************/
int bit_reverse(int k, int r)
{
    int l = 0;                                 // Accumulate the results;
    for(int i=0; i<r; i++) {                   // Loop on every bit;
        l = (l << 1) + (k & 1);                // Test less signficant bit and add;
        k >>= 1;                               // Test next bit;
    }
    return l;
}


/**************************************************************************************************
 * Function: iterative_fft
 *   Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
 *   O(N log_2(N)) complexity, and since there are less function calls, it will probably be
 *   marginally faster than the recursive versions.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector of
 *     a power of two length, or it will fail. No checks on this are made.
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void iterative_fft(Complex x[], Complex X[], int N)
{
    Complex W, Wkn, w;                         // Twiddle factors;
    int r, l, p, q, step;

    r = (int) floor(log(N)/log(2));            // Number of bits;
    for(int k=0; k<N; k++) {
        l = bit_reverse(k, r);                 // Reorder the vector according to the
        X[l] = x[k];                           //   bit-reversed order;
    }

    step = 1;                                  // Auxiliary for computation of twiddle factors;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            W = cexpn(-M_PI/step);             // Twiddle factors;
            Wkn = cmplx(1, 0);
            for(int n=0; n<step; n++) {
                p = l + n;
                q = p + step;
                w = cmul(X[q], Wkn);           // Recombine results;
                X[q] = csub(X[p], w);
                w = cmul(X[p], cmplx(2, 0));
                X[p] = csub(w, X[q]);
                Wkn = cmul(Wkn, W);            // Update twiddle factors;
            }
        }
        step <<= 1;
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, itime;
    int n;

    // Start by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |\n");
    printf("+---------+---------+---------+---------+---------+---------+\n");

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int r=5; r<11; r++) {

        // Compute the average execution time:
        n = exp2(r);
        dtime = time_it(direct_ft, n, REPEAT);
        rtime = time_it(recursive_fft, n, REPEAT);
        itime = time_it(iterative_fft, n, REPEAT);

        // Print the results:
        printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, r*n, dtime, rtime, itime);
    }
    printf("+---------+---------+---------+---------+---------+---------+\n");

    return 0;
}
