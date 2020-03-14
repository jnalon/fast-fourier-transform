/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program was compiled and tested in Linux using Mono command line tools. It might need some
 * adaptation to be compiled with Visual Studio (eg. creating a project), but it can be compiles as
 * is with Windows command line tools:
 *
 * $ mcs anyfft.cs
 *
 * This will generate a file called 'anyfft.exe', that can be run with the command (remember to
 * change permission to execute):
 *
 * $ ./anyfft.exe
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
     * Auxiliary Method: time_it
     *   This function calls a Fast Fourier Transform function repeatedly a certain number of
     *   times, measure execution time and average it.
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
     **********************************************************************************************/
    public static double TimeIt(DFT f, int size, int repeat)
    {
        Complex[] x = new Complex[size];
        Complex[] X;

        for(int i=0; i<size; i++)
            x[i] = new Complex(i, 0);
        Stopwatch dsw = Stopwatch.StartNew();
        for(int j=0; j<repeat; j++)
            X = f(x);
        dsw.Stop();
        double dtime = ((double) dsw.ElapsedMilliseconds) / ((double)(1000*repeat));
        return dtime;
    }

    /**********************************************************************************************
     * Method: DirectFT
     *   Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that
     *   has O(N^2) complexity.
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
        int N = x.Length;
        Complex[] X = new Complex[N];
        Complex W = Complex.exp(-2*Math.PI/N);         // Initializes twiddle factors;
        Complex Wk = new Complex(1, 0);

        for(int k=0; k<N; k++) {
            X[k] = new Complex();                      // Accumulates the results;
            Complex Wkn = new Complex(1, 0);
            for(int n=0; n<N; n++) {
                X[k] = X[k] + Wkn*x[n];
                Wkn = Wkn * Wk;                        // Update twiddle factor;
            }
            Wk = Wk * W;
        }
        return X;
    }


    /**********************************************************************************************
     * Method: Factor
     *   This method finds the smallest prime factor of a given number. If the argument is prime
     *   itself, then it is the return value.
     *
     * Parameters:
     *   n
     *     Number to be inspected.
     *
     * Returns:
     *   The smallest prime factor, or the number itself if it is already a prime.
     **********************************************************************************************/
    public static int Factor(int n)
    {
        int rn = n/2;                              // Search up to half the number;
        for(int i=2; i<=rn; i++)
            if (n%i == 0) return i;                // If remainder is zero, a factor is found;
        return n;
    }



    /**********************************************************************************************
     * Method: RecursiveFFT
     *   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This
     *   has smallest complexity than the direct FT, though the exact value is difficult to
     *   compute.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. It must be a composite number, or else the
     *     computation will be defered to the direct FT, and there will be no efficiency gain.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] RecursiveFFT(Complex[] x)
    {
        int N = x.Length;
        Complex[] X = new Complex[N];
        for(int n=0; n<N; n++)
            X[n] = new Complex();

        int N1 = Factor(N);                            // Smallest prime factor of length;
        if (N1==N)                                     // If the length is prime itself,
            return DirectFT(x);                        //   transform is given by the direct form;
        else {
            int N2 = N / N1;                           // Decompose in two factors, N1 being prime;

            Complex[] xj = new Complex[N2];            // Allocates memory for subsequences;

            Complex W = Complex.exp(-2*Math.PI/N);     // Twiddle factor;
            Complex Wj = new Complex(1, 0);
            for(int j=0; j<N1; j++) {                  // Computes every subsequence of size N2;
                for(int n=0; n<N2; n++)
                    xj[n] = x[n*N1+j];                 // Creates the subsequence;
                Complex[] Xj = RecursiveFFT(xj);       // Computes the DFT of the subsequence;
                Complex Wkj = new Complex(1, 0);
                for(int k=0; k<N; k++) {               // Recombine results;
                    X[k] = X[k] + Xj[k%N2] * Wkj;
                    Wkj = Wkj * Wj;                    // Update twiddle factors;
                }
                Wj = Wj * W;
            }
            return X;
        }
    }


    /**********************************************************************************************
     * Main Method:
     **********************************************************************************************/
    public static void Main()
    {
        int[] SIZES = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
        int REPEAT = 500;                      // Number of executions to compute average time;

        // Starts by printing the table with time comparisons:
        Console.WriteLine("+---------+---------+---------+---------+");
        Console.WriteLine("|    N    |   N^2   | Direta  | Recurs. |");
        Console.WriteLine("+---------+---------+---------+---------+");

        // Try it with vectors with the given sizes:
        for(int i=0; i<SIZES.Length; i++) {

            // Computes the average execution time:
            int n = SIZES[i];
            double dtime = TimeIt(DirectFT, n, REPEAT);
            double rtime = TimeIt(RecursiveFFT, n, REPEAT);

            // Print the results:
            string results = String.Format("| {0,7} | {1,7} | {2,7:F4} | {3,7:F4} |",
                n, n*n, dtime, rtime);
            Console.WriteLine(results);
        }

        Console.WriteLine("+---------+---------+---------+---------+");
    }

}
