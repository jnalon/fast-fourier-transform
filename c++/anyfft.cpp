/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C++ compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ g++ -o anyfft anyfft.cpp -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./anyfft
 *
 * Obs.: Technically, the power of C++ resides in the object orientation facilities. This program,
 *   however, doesn't use a lot of it, given its simplicity: it mainly operates over a vector. In a
 *   object orientation environment (a big project, for instance), maybe the best way to do it was
 *   to create a complex vector class and make the FFT a method of it. The same could be said of a
 *   number of other resources such as arrays and libraries, but we'll keep it simple here.
 **************************************************************************************************/

/**************************************************************************************************
 Include necessary libraries:
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
 Small class to operate with complex numbers:
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
    Complex W = cexpn(-2*M_PI/N);              // Initialize twiddle factors;
    Complex Wk = Complex(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = Complex();                      // Accumulate the results;
        Complex Wkn = Complex(1, 0);           // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            X[k] = X[k] + Wkn*x[n];
            Wkn = Wkn * Wk;                    // Update twiddle factor;
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
    int rn = n/2;                              // Search up to half the number;
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
 *     The vector of which the FFT will be computed. It must be a composite number, or else the
 *     computation will be defered to the direct FT, and there will be no efficiency gain.
 *   X
 *     The vector that will receive the results of the computation. It needs to be allocated prior
 *     to the function call;
 *   N
 *     The number of elements in the vector.
 **************************************************************************************************/
void recursive_fft(Complex x[], Complex X[], int N)
{
    int N1 = factor(N);                        // Smallest prime factor of length;
    if(N1==N)                                  // If the length is prime itself,
        direct_ft(x, X, N);                    //   the transform is given by the direct form;
    else {
        int N2 = N / N1;                       // Decompose in two factors, N1 being prime;

        Complex *xj = new Complex[N2];         // Allocate memory for subsequences;
        Complex *Xj = new Complex[N2];         //   and their transforms;

        Complex W = cexpn(-2*M_PI/N);          // Twiddle factor;
        Complex Wj = Complex(1, 0);
        for(int j=0; j<N1; j++) {              // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++) {
                xj[n] = x[n*N1+j];             // Create the subsequence;
                Xj[n] = Complex(0, 0);         // Initialize the transform;
            }
            recursive_fft(xj, Xj, N2);         // Compute the DFT of the subsequence;
            Complex Wkj = Complex(1, 0);
            for(int k=0; k<N; k++) {
                X[k] = X[k] + Xj[k%N2] * Wkj;  // Recombine results;
                Wkj = Wkj * Wj;                // Update twiddle factors;
            }
            Wj = Wj * W;
        }

        delete Xj;                              // Clean-up used memory;
        delete xj;
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    int SIZES[] = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };

    // Start by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | Direct  | Recurs. |" << endl;
    cout << "+---------+---------+---------+---------+" << endl;

    // Try it with vectors with the given sizes:
    for(int i=0; i<8; i++) {

        // Compute the average execution time:
        int n = SIZES[i];
        float dtime = time_it(direct_ft, n, REPEAT);
        float rtime = time_it(recursive_fft, n, REPEAT);

        // Print the results:
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) << setprecision(7) << dtime << " ";
        cout << "| " << setw(7) << setprecision(7) << rtime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------+" << endl;
    return 0;
}
