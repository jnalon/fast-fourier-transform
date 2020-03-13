/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ g++ -o fft fft.cpp -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fft
 *
 * Obs.: Technically, the power of C++ resides in the object orientation facilities. This program,
 *   however, doesn't use a lot of it, given its simplicity: it mainly operates over a vector. In a
 *   object orientation environment (a big project, for instance), maybe the best way to do it was
 *   to create a complex vector class and make the FFT a method of it. The same could be said of a
 *   number of other resources such as arrays and libraries, but we'll keep it simple here.
 **************************************************************************************************/

/**************************************************************************************************
 Includes necessary libraries:
 **************************************************************************************************/
#include <iostream>                            // Input and Output;
#include <iomanip>                             // I/O Manipulation;
#include <array>                               // Deals with arrays;
#include <cmath>                               // Math Functions;
#include <chrono>                              // Time measurement;

using namespace std;


/**************************************************************************************************
 Definitions:
 **************************************************************************************************/
#define REPEAT 500                             // Number of executions to compute average time;


/**************************************************************************************************
 Small class to operate with complex numbers
 **************************************************************************************************/
class Complex {
    public:
        float r;                               // Real part;
        float i;                               // Imaginary part;
        Complex();                             // Constructors;
        Complex(float re, float im);
        void set(float re, float im);
        void set(Complex c);
        Complex operator+(Complex c);          // Addition (overload + operator);
        Complex operator-(Complex c);          // Subtraction (overload - operator);
        Complex operator*(Complex c);          // Product (overload * operator);
        Complex operator*(float a);            // Product with a scalar;
        Complex cexp();                        // Complex exponential;
};

Complex::Complex() {                           // Constructor;
    r = 0.0;                                   // Real part;
    i = 0.0;                                   // Imaginary part;
}

Complex::Complex(float re, float im) {         // Constructor;
    r = re;                                    // Real part;
    i = im;                                    // Imaginary part;
}

void Complex::set(float re, float im) {
    r = re;                                    // Real part;
    i = im;                                    // Imaginary part;
}

void Complex::set(Complex c) {
    r = c.r;                                   // Real part;
    i = c.i;                                   // Imaginary part;
}

Complex Complex::operator+(Complex c) {
    return Complex(r + c.r, i + c.i);
}

Complex Complex::operator-(Complex c) {
    return Complex(r - c.r, i - c.i);
}

Complex Complex::operator*(Complex c) {
    return Complex(r*c.r - i*c.i, r*c.i + i*c.r);
}

Complex Complex::operator*(float a) {
    return Complex(a*r, a*i);
}

Complex Complex::cexp() {
    return Complex(exp(r)*cos(i), exp(r)*sin(i));
}

Complex cexpn(float a) {                       // Convenience function to compute the exponential;
    return Complex(cos(a), sin(a));
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

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = Complex();
    auto t0 = chrono::steady_clock::now();     // Start counting time;
    for(int j=0; j<repeat; j++)
        (*f)(x, X, size);
    auto t1 = chrono::steady_clock::now();     // End of time measuring;
    return chrono::duration_cast<chrono::seconds>(t1 - t0).count() / (float) repeat;
}


/**************************************************************************************************
 * Function: direct_ft
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
    Complex W = cexpn(-2*M_PI/N);              // Initializes twiddle factors;
    Complex Wk = Complex(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = Complex();                      // Accumulates the results;
        Complex Wkn = Complex(1, 0);           // Initializes twiddle factors;
        for(int n=0; n<N; n++) {
            X[k] = X[k] + Wkn*x[n];
            Wkn = Wkn * Wk;                    // Update twiddle factor;
        }
        Wk = Wk * W;
    }
    return;
}


/**************************************************************************************************
 * Function: recursive_fft
 *   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This has
 *   O(N log_2(N)) complexity.
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
    if(N==1) {                                 // A length-1 vector is its own FT;
        X[0] = x[0];
        return;
    }
    int N2 = N >> 1;
    Complex *xe = new Complex[N2];             // Allocates memory for computation;
    Complex *xo = new Complex[N2];
    Complex *Xe = new Complex[N2];
    Complex *Xo = new Complex[N2];
    for(int k=0; k<N2; k++) {                  // Split even and odd samples;
        xe[k] = x[k<<1];
        xo[k] = x[(k<<1)+1];
    }
    recursive_fft(xe, Xe, N2);                 // Transform of even samples;
    recursive_fft(xo, Xo, N2);                 // Transform of odd samples;
    Complex W = cexpn(-2*M_PI/N);              // Twiddle factors;
    Complex Wk = Complex(1, 0);
    for(int k=0; k<N2; k++) {
        Complex w = Wk * Xo[k];                // Recombine results;
        X[k] = Xe[k] + w;
        X[k+N2] = Xe[k] - w;
        Wk = Wk * W;                           // Update twiddle factors;
    }
    delete Xo;                                 // Releasing memory of intermediate vectors;
    delete Xe;
    delete xo;
    delete xe;
    return;
}


/**************************************************************************************************
 * Function: bit_reverse
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
 * Function: iterative_fft
 *   Computes the Fast Fourier Ttransform using an iterative in-place decimation in time algorithm.
 *   This has O(N log_2(N)) complexity, and since there are less function calls, it will probably
 *   be marginally faster than the recursive versions.
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
    int r = (int) floor(log(N)/log(2));
    for(int k=0; k<N; k++) {
        int l = bit_reverse(k, r);
        X[l] = x[k];
    }
    int step = 1;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            Complex W = cexpn(-M_PI/step);
            Complex Wkn = Complex(1, 0);
            for(int n=0; n<step; n++) {
                int p = l + n;
                int q = p + step;
                X[q] = X[p] - Wkn * X[q];
                X[p] = X[p]*2 - X[q];
                Wkn = Wkn * W;
            }
        }
        step <<= 1;
    }
    return;
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    // Starts by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |" << endl;
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int r=5; r<11; r++) {

        // Computes the average execution time:
        int n = (int) exp2(r);
        float dtime = time_it(direct_ft, n, REPEAT);
        float rtime = time_it(recursive_fft, n, REPEAT);
        float itime = time_it(iterative_fft, n, REPEAT);

        // Print the results:
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) <<   r*n << " ";
        cout << "| " << setw(7) << setprecision(7) << dtime << " ";
        cout << "| " << setw(7) << setprecision(7) << rtime << " ";
        cout << "| " << setw(7) << setprecision(7) << itime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    return 0;
}
