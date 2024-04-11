/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This file implements methods for testing and time measuring.
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;
using System.Diagnostics;                      // Time measuring;
// using Complex;


class Test
{
    /// <summary>
    /// Initializes a complex vector. This method creates and initializes a vector of complex
    /// numbers with bogus data for testing purposes only.
    /// <param name="size">Number of elements in the vector.</param>
    /// <returns> The initialized vector.</returns>
    /// </summary>
    public static Complex[] InitializeVector(int size)
    {
        Complex[] x = new Complex[size];
        for(int i=0; i<size; i++)                      // Initialize the vector;
            x[i] = new Complex(i, 0);
        return x;
    }


    /// <summary>
    /// Pretty printing of an array of complex numbers, used to inspect results.
    /// <param name="f">
    ///     Function to be called, with the given prototype. The parameter is a complex vector to
    ///     be transformed, and the return value is a complex vector with the coefficients of the
    ///     transform;
    /// </param>
    /// <param name="size">
    ///     Number of elements in the vector on which the transform will be applied;
    /// </param>
    /// </summary>
    public static void TestIt(FFT.DFT f, int size)
    {
        Complex[] x = InitializeVector(size);
        Complex[] X = f(x);
        Console.WriteLine(String.Format("N = {0,2} | Input | Output", size));
        for(int i=0; i<size; i++)
             Console.WriteLine(String.Format("  {0, 2}: {1} | {2}", i, x[i], X[i]));
        Console.WriteLine("------------------------------");
    }


    /// <summary>
    /// Pretty printing of an array of complex numbers, used to inspect results.
    /// <param name="f">
    ///     Function to be called, with the given prototype. The parameter is a complex vector to
    ///     be transformed, and the return value is a complex vector with the coefficients of the
    ///     transform;
    /// </param>
    /// <param name="size">
    ///     Number of elements in the vector on which the transform will be applied;
    /// </param>
    /// <param name="repeat">Number of times the function will be called.</param>
    /// <returns>
    ///     The average execution time for that function with a vector of the given size.
    /// </returns>
    /// </summary>
    public static double TimeIt(FFT.DFT f, int size, int repeat)
    {
        Complex[] x = InitializeVector(size);
        Stopwatch dsw = Stopwatch.StartNew();          // Start a timer;
        for(int j=0; j<repeat; j++)                    // Repeated calls;
            f(x);
        dsw.Stop();
        return ((double) dsw.ElapsedMilliseconds) / ((double)(1000*repeat));
    }

}
