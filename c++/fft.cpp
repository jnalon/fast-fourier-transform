/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C++ compiler, just remember to link the math library. In my box, I used the command:
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
 Include necessary libraries:
 **************************************************************************************************/
#include <iostream>                            // Input and Output;
#include <iomanip>                             // I/O Manipulation;
#include <vector>                              // Deals with arrays;
#include <cmath>                               // Math Functions;
#include <chrono>                              // Time measurement;

using namespace std;


/**************************************************************************************************
 Definitions:
 **************************************************************************************************/
const int repeats = 500;                       // Number of executions to compute average time;


/**************************************************************************************************
 Small class to operate with complex numbers:
 **************************************************************************************************/
class Complex {
    public:
        float r;                               // Real part;
        float i;                               // Imaginary part;

    Complex() {                                // Constructor;
        r = 0.0;
        i = 0.0;
    }

    Complex(float re, float im) {              // Constructor;
        r = re;
        i = im;
    }

    Complex operator+(Complex c) {
        return Complex(r + c.r, i + c.i);
    }

    Complex operator-(Complex c) {
        return Complex(r - c.r, i - c.i);
    }

    Complex operator*(Complex c) {
        return Complex(r*c.r - i*c.i, r*c.i + i*c.r);
    }

    Complex operator*(float a) {
        return Complex(a*r, a*i);
    }
};

Complex cexpn(float a) {                       // Convenience function to compute the exponential;
    return Complex(cos(a), sin(a));
}

std::ostream& operator <<( std::ostream& out_stream, const Complex x)
{
    return out_stream << "(" << x.r << ", " << x.i << ")";
}


/**************************************************************************************************
 * Auxiliary function: complex_show
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 **************************************************************************************************/
void complex_show(vector<Complex> x)
{
    for (int i=0; i<x.size(); i++)
        std::cout << x[i] << std::endl;
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
float time_it(vector<Complex> (*f)(vector<Complex>), int size, int repeat)
{
    vector<Complex> x(size), X;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = Complex(j, 0);
    auto t0 = chrono::steady_clock::now();     // Start counting time;
    for(int j=0; j<repeat; j++)
        X = (*f)(x);
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
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
vector<Complex> direct_ft(vector<Complex> x)
{
    int N = x.size();
    vector<Complex> X(N);
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
    return X;
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
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
vector<Complex> recursive_fft(vector<Complex> x)
{
    int N = x.size();
    if(N == 1) {                                       // A length-1 vector is its own FT;
        vector<Complex> X = x;
        return X;
    } else {
        int N2 = N >> 1;
        vector<Complex> xe(N2);                        // Allocate memory for computation;
        vector<Complex> xo(N2);
        vector<Complex> X(N);

        for(int k=0; k<N2; k++) {                      // Split even and odd samples;
            xe[k] = x[k<<1];
            xo[k] = x[(k<<1)+1];
        }
        vector<Complex> Xe = recursive_fft(xe);        // Transform of even samples;
        vector<Complex> Xo = recursive_fft(xo);        // Transform of odd samples;

        Complex W = cexpn(-2*M_PI/N);                  // Twiddle factors;
        Complex Wk = Complex(1, 0);
        for(int k=0; k<N2; k++) {
            Complex w = Wk * Xo[k];                    // Recombine results;
            X[k] = Xe[k] + w;
            X[k+N2] = Xe[k] - w;
            Wk = Wk * W;                               // Update twiddle factors;
        }
        return X;
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
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
vector<Complex> iterative_fft(vector<Complex> x)
{
    int N = x.size();
    int r = (int) floor(log2(N));              // Number of bits;
    vector<Complex> X(N);
    for(int k=0; k<N; k++) {
        int l = bit_reverse(k, r);             // Reorder the vector according to the
        X[l] = x[k];                           //   bit-reversed order;
    }

    int step = 1;                              // Auxiliary for computation of twiddle factors;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            Complex W = cexpn(-M_PI/step);     // Twiddle factors;
            Complex Wkn = Complex(1, 0);
            for(int n=0; n<step; n++) {
                int p = l + n;
                int q = p + step;
                X[q] = X[p] - Wkn * X[q];      // Recombine results;
                X[p] = X[p]*2 - X[q];
                Wkn = Wkn * W;                 // Update twiddle factors;
            }
        }
        step <<= 1;
    }
    return X;
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    // Start by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  |" << endl;
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int r=5; r<11; r++) {

        // Compute the average execution time:
        int n = (int) exp2(r);
        float dtime = time_it(direct_ft, n, repeats);
        float rtime = time_it(recursive_fft, n, repeats);
        float itime = time_it(iterative_fft, n, repeats);

        // Print the results:
        cout << fixed;
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) <<   r*n << " ";
        cout << "| " << setw(7) << setprecision(4) << dtime << " ";
        cout << "| " << setw(7) << setprecision(4) << rtime << " ";
        cout << "| " << setw(7) << setprecision(4) << itime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    return 0;
}
