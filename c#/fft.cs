/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program was compiled and tested in Linux using Mono command line tools. It might need some
 * adaptation to be compiled with Visual Studio (eg. creating a project), but it can be compiles as
 * is with Windows command line tools:
 *
 * $ mcs fft.cs
 *
 * This will generate a file called 'fft.exe', that can be run with the command (remember to change
 * permission to execute):
 *
 * $ ./fft.exe
 **************************************************************************************************/

/**************************************************************************************************
 Include necessary libraries:
 **************************************************************************************************/
using System;                                  // Input and output and standard library;
using System.Diagnostics;                      // Time measuring;


/**************************************************************************************************
 Small class to operate with complex numbers:
 **************************************************************************************************/
public class Complex
{
    public double re;                          // Real part;
    public double im;                          // Imaginary part;

    public Complex()
    {
        this.re = 0;
        this.im = 0;
    }

    public Complex(double real, double imag)
    {
        this.re = real;
        this.im = imag;
    }

    // Sum of two complex numbers:
    public static Complex operator +(Complex z1, Complex z2)
    {
        return new Complex(z1.re + z2.re, z1.im + z2.im);
    }

    // Difference of two complex numbers:
    public static Complex operator -(Complex z1, Complex z2)
    {
        return new Complex(z1.re - z2.re, z1.im - z2.im);
    }

    // Product of two complex numbers:
    public static Complex operator *(Complex z1, Complex z2)
    {
        return new Complex(z1.re*z2.re - z1.im*z2.im, z1.re*z2.im + z1.im*z2.re);
    }

    // Product with a scalar:
    public static Complex operator *(Complex z, double x)
    {
        return new Complex(z.re*x, z.im*x);
    }

    // Complex Exponential:
    public static Complex exp(double a)
    {
        return new Complex(Math.Cos(a), Math.Sin(a));
    }

    // String representation for printing:
    public override string ToString()
    {
        return(String.Format("{0} + {1}i", re, im));
    }
}


/**************************************************************************************************
 Main program and routines:
 **************************************************************************************************/
class FastFourierTransform
{
    // Prototype to the DFT function, to be used at the timing function:
    public delegate Complex[] DFT(Complex[] x);


    /**********************************************************************************************
     * Auxiliary Method: ComplexShow
     *   Pretty printing of an array of complex numbers, used to inspect results.
     *
     * Parameters:
     *   x
     *     A vector of complex numbers, according to the definition above;
     **********************************************************************************************/
    public static void ComplexShow(Complex[] x)
    {
        for(int k=0; k<x.Length; k++)
            Console.WriteLine(x[k].ToString());
    }


    /**********************************************************************************************
     * Auxiliary Method: TimeIt
     *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
     *
     * Parameters:
     *  f
     *    Function to be called, with the given prototype. The parameter is a complex vector to be
     *    transformed, and the return value is a complex vector with the coefficients of the
     *    transform;
     *  size
     *    Number of elements in the vector on which the transform will be applied;
     *  repeat
     *    Number of times the function will be called.
     *
     * Returns:
     *   The average execution time for that function with a vector of the given size.
     **********************************************************************************************/
    public static double TimeIt(DFT f, int size, int repeat)
    {
        Complex[] x = new Complex[size];

        for(int i=0; i<size; i++)                      // Initialize the vector;
            x[i] = new Complex(i, 0);
        Stopwatch dsw = Stopwatch.StartNew();          // Start a timer;
        for(int j=0; j<repeat; j++)                    // Repeated calls;
            f(x);
        dsw.Stop();
        return ((double) dsw.ElapsedMilliseconds) / ((double)(1000*repeat));
    }


    /**********************************************************************************************
     * Method: DirectFT
     *   Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
     *   complexity.
     *
     * Parameters:
     *   x
     *     The vector of which the DFT will be computed. Given the nature of the implementation,
     *     there is no restriction on the size of the vector, although it will almost always be
     *     called with a power of two size to give a fair comparison;
     *
     * Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] DirectFT(Complex[] x)
    {
        int N = x.Length;                              // Length of the vector;
        Complex[] X = new Complex[N];                  // Accumulate the results;
        Complex W = Complex.exp(-2*Math.PI/N);         // Initialize twiddle factors;
        Complex Wk = new Complex(1, 0);

        for(int k=0; k<N; k++) {                       // Compute the kth coefficient;
            X[k] = new Complex();                      // Accumulate the results;
            Complex Wkn = new Complex(1, 0);
            for(int n=0; n<N; n++) {                   //   Operate the summation;
                X[k] = X[k] + Wkn*x[n];                //     Compute every term;
                Wkn = Wkn * Wk;                        // Update twiddle factor;
            }
            Wk = Wk * W;
        }
        return X;
    }


    /**********************************************************************************************
     * Method: RecursiveFFT
     *   Fast Fourier Transform using a recursive decimation in time algorithm. This has
     *   O(N log_2(N)) complexity.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. This should always be called with a vector
     *     of a power of two length, or it will fail. No checks on this are made.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] RecursiveFFT(Complex[] x)
    {
        int N = x.Length;

        if (N==1)                                      // A length-1 vector is its own FT;
            return x;
        else {
            int N2 = N >> 1;

            Complex[] xe = new Complex[N2];            // Allocate memory for computation;
            Complex[] xo = new Complex[N2];
            Complex[] X = new Complex[N];

            for(int k=0; k<N2; k++) {                  // Split even and odd samples;
                xe[k] = x[k<<1];
                xo[k] = x[(k<<1)+1];
            }
            Complex[] Xe = RecursiveFFT(xe);           // Transform of even samples;
            Complex[] Xo = RecursiveFFT(xo);           // Transform of odd samples;

            Complex W = Complex.exp(-2*Math.PI/N);     // Twiddle factors;
            Complex Wk = new Complex(1, 0);
            for(int k=0; k<N2; k++) {
                Complex w = Wk * Xo[k];                // Recombine results;
                X[k] = Xe[k] + w;
                X[k+N2] = Xe[k] - w;
                Wk = Wk * W;                           // Update twiddle factors;
            }
            return X;
        }
    }


    /**********************************************************************************************
     * Method: BitReverse
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
     **********************************************************************************************/
    private static int BitReverse(int k, int r)
    {
        int l = 0;                                     // Accumulate the results;
        for(int i=0; i<r; i++) {                       // Loop on every bit;
            l = (l << 1) + (k & 1);                    // Test less signficant bit and add;
            k = (k >> 1);                              // Test next bit;
        }
        return l;
    }


    /**********************************************************************************************
     * Function: IterativeFFT
     *   Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
     *   O(N log_2(N)) complexity, and since there are less function calls, it will probably be
     *   marginally faster than the recursive versions.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. This should always be called with a vector
     *     of a power of two length, or it will fail. No checks on this are made.
     *
     * Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] IterativeFFT(Complex[] x)
    {
        int N = x.Length;
        Complex[] X = new Complex[N];

        int r = (int) Math.Round(Math.Log(N, 2));              // Number of bits;
        for(int k=0; k<N; k++) {
            int l = BitReverse(k, r);                          // Reorder the vector according to
            X[l] = x[k];                                       //   the bit-reversed order;
        }

        int step = 1;                                          // Computation of twiddle factors;
        for(int k=0; k<r; k++) {
            for(int l=0; l<N; l+=2*step) {
                Complex W = Complex.exp(-Math.PI/step);        // Twiddle factors;
                Complex Wkn = new Complex(1, 0);
                for(int n=0; n<step; n++) {
                    int p = l + n;
                    int q = p + step;
                    X[q] = X[p] - Wkn * X[q];                  // Recombine results;
                    X[p] = X[p]*2 - X[q];
                    Wkn = Wkn * W;                             // Update twiddle factors;
                }
            }
            step <<= 1;
        }

        return X;
    }


    /**********************************************************************************************
     * Main Method:
     **********************************************************************************************/
    public static void Main()
    {
        int REPEAT = 500;                      // Number of executions to compute average time;

        // Start by printing the table with time comparisons:
        Console.WriteLine("+---------+---------+---------+---------+---------+---------+");
        Console.WriteLine("|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |");
        Console.WriteLine("+---------+---------+---------+---------+---------+---------+");

        // Try it with vectors with size ranging from 32 to 1024 samples:
        for(int r=5; r<11; r++) {

            // Compute the average execution time:
            int n = (int) Math.Pow(2, r);
            double dtime = TimeIt(DirectFT, n, REPEAT);
            double rtime = TimeIt(RecursiveFFT, n, REPEAT);
            double itime = TimeIt(IterativeFFT, n, REPEAT);

            // Print the results:
            string results = String.Format("| {0,7} | {1,7} | {2,7} | {3,7:F4} | {4,7:F4} | {5,7:F4} |",
                n, n*n, n*r, dtime, rtime, itime);
            Console.WriteLine(results);
        }

        Console.WriteLine("+---------+---------+---------+---------+---------+---------+");
    }

}
