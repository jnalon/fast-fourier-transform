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
 Includes necessary libraries:
 **************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


/**************************************************************************************************
 Data structure to deal with complex numbers. There is, of course, an external library that
 implements operations for comple numbers, but I'm defining this here to reduce dependencies.
 **************************************************************************************************/
typedef struct {
    float r, i;                                // Real and Imaginary parts of a complex number;
} Complex;

#define REPEAT 500                             // Number of executions to compute average time;


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
    int i;

    for (i=0; i<n; ++i)
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
 *  n
 *    Number of elements in the vector on which the transform will be applied;
 *  repeat
 *    Number of times the function will be called.
 *
 * Returns:
 *   The average execution time for that function with a vector of the given size.
 **************************************************************************************************/
float time_it(void (*f)(Complex *, Complex *, int), int n, int repeat)
{
    Complex x[1024], X[1024];                  // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;
    int j;

    for(j=0; j<n; ++j) {                       // Initialize the vector;
        x[j].r = j;
        x[j].i = 0;
    }
    t0 = clock();                              // Start counting time;
    for(j=0; j<repeat; ++j)
        (*f)(x, X, n);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeat;
    return t1;
}


/**************************************************************************************************
 * function: direct_ft
 *   Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that has
 *   O(N^2) complexity.
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
    int k, n;

    W = cexpn(-2*M_PI/N);                      // Initializes twiddle factors;
    Wk = cmplx(1, 0);
    for(k=0; k<N; ++k) {
        X[k] = cmplx(0, 0);                    // Accumulates the results;
        Wkn = cmplx(1, 0);                     // Initializes twiddle factors;
        for(n=0; n<N; ++n) {
            w = cmul(x[n], Wkn);
            X[k] = cadd(X[k], w);
            Wkn = cmul(Wkn, Wk);               // Updates twiddle factors;
        }
        Wk = cmul(Wk, W);
    }
    return;
}


/**************************************************************************************************
 * function: recursive_fft
 *   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This has
 *   O(N log_2(N)) complexity. This implementation uses native Python lists.
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
    int k, N2;

    if(N==1) {                                 // A length-1 vector is its own FT;
        X[0] = x[0];
        return;
    }
    N2 = N >> 1;
    xe = (Complex *) malloc(sizeof(Complex)*N2);       // Allocates memory for computation;
    xo = (Complex *) malloc(sizeof(Complex)*N2);
    Xe = (Complex *) malloc(sizeof(Complex)*N2);
    Xo = (Complex *) malloc(sizeof(Complex)*N2);
    for(k=0; k<N2; ++k) {                              // Separates even and odd samples;
        xe[k] = x[k<<1];
        xo[k] = x[(k<<1)+1];
    }
    recursive_fft(xe, Xe, N2);                         // Transform of even samples;
    recursive_fft(xo, Xo, N2);                         // Transform of odd samples;
    W = cexpn(-2*M_PI/N);                              // Twiddle factors;
    Wk = cmplx(1, 0);
    for(k=0; k<N2; ++k) {
        w = cmul(Xo[k], Wk);
        X[k] = cadd(Xe[k], w);                         // Recombine results;
        X[k+N2] = csub(Xe[k], w);
        Wk = cmul(Wk, W);                              // Update twiddle factors;
    }
    free(Xo);                                          // Releasing memory of intermediate vectors;
    free(Xe);
    free(xo);
    free(xe);
    return;
}


/**************************************************************************************************
 * function: bit_reverse
 *   Computes the bit-reversed function of an integer number.
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
    int l, i;

    l = 0;                                     // Accumulates the results;
    for(i=0; i<r; i++) {                       // Loop on every bit;
        l = (l << 1) + (k & 1);                // Tests less signficant bit and add;
        k >>= 1;                               // Tests next bit;
    }
    return l;
}


/**************************************************************************************************
 * function: interactive_fft
 *   Computes the Fast Fourier Ttransform using an interactive in-place decimation in time
 *   algorithm. This has O(N log_2(N)) complexity, and since there are less function calls, it will
 *   probably be marginally faster than the recursive versions.
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
void interactive_fft(Complex x[], Complex X[], int N)
{
    Complex W, Wkn, w;                         // Twiddle factors;
    int r, k, l, n, p, q, step;

    r = (int) floor(log(N)/log(2));            // Number of bits;
    for(k=0; k<N; ++k) {                       // Reorder the vector according to the
        l = bit_reverse(k, r);                 //   bit-reversed order;
        X[l] = x[k];
    }
    step = 1;                                  // Auxiliary for computation of twiddle factors;
    for(k=0; k<r; ++k) {
        for(l=0; l<N; l+=2*step) {
            W = cexpn(-M_PI/step);             // Twiddle factors;
            Wkn = cmplx(1, 0);
            for(n=0; n<step; ++n) {
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
    return;
}


/**************************************************************************************************
 Funcao Principal.
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, itime;
    int r, n;

    // Starts by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |\n");
    printf("+---------+---------+---------+---------+---------+---------+\n");

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(r=5; r<11; ++r) {

        // Computes the average execution time:
        n = exp2(r);
        dtime = time_it(direct_ft, n, REPEAT);
        rtime = time_it(recursive_fft, n, REPEAT);
        itime = time_it(interactive_fft, n, REPEAT);

        // Print the results:
        printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, r*n, dtime, rtime, itime);
    }
    printf("+---------+---------+---------+---------+---------+---------+\n");

    return 0;
}
