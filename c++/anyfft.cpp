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
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
vector<Complex> recursive_fft(vector<Complex> x)
{
    int N = x.size();
    int N1 = factor(N);                        // Smallest prime factor of length;
    if(N1 == N) {                              // If the length is prime itself,
        vector<Complex> X = direct_ft(x);      //   the transform is given by the direct form;
        return X;
    } else {
        int N2 = N / N1;                       // Decompose in two factors, N1 being prime;
        vector<Complex> xj(N2);                // Allocate memory for subsequences;
        vector<Complex> X(N);

        Complex W = cexpn(-2*M_PI/N);          // Twiddle factor;
        Complex Wj = Complex(1, 0);
        for(int j=0; j<N1; j++) {              // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++)
                xj[n] = x[n*N1+j];             // Create the subsequence;
            vector<Complex> Xj = recursive_fft(xj);         // Compute the DFT of the subsequence;
            Complex Wkj = Complex(1, 0);
            for(int k=0; k<N; k++) {
                X[k] = X[k] + Xj[k%N2] * Wkj;  // Recombine results;
                Wkj = Wkj * Wj;                // Update twiddle factors;
            }
            Wj = Wj * W;
        }
        return X;
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    vector<int> sizes = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };

    // Start by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | Direct  | Recurs. |" << endl;
    cout << "+---------+---------+---------+---------+" << endl;

    // Try it with vectors with the given sizes:
    for(int& n : sizes) {

        // Compute the average execution time:
        float dtime = time_it(direct_ft, n, repeats);
        float rtime = time_it(recursive_fft, n, repeats);

        // Print the results:
        cout << fixed;
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) << setprecision(4) << dtime << " ";
        cout << "| " << setw(7) << setprecision(4) << rtime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------+" << endl;
    return 0;
}
