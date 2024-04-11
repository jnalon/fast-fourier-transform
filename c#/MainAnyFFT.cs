/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 **************************************************************************************************
 * This program was compiled and tested in Linux using Mono command line tools. It might need some
 * adaptation to be compiled with Visual Studio (eg. creating a project), but it can be compiles as
 * is with Windows command line tools:
 *
 * $ mcs MainAnyFFT.cs Complex.cs FFT.cs Test.cs
 *
 * This will generate a file called 'MainAnyFFT.exe', that can be run with the command:
 *
 * $ ./MainAnyFFT.exe
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;


class MainFFT
{
    public static void Main()
    {
        int[] SIZES = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
        int REPEAT = 500;                      // Number of executions to compute average time;

        // Starts by printing the table with time comparisons:
        Console.WriteLine("+---------+---------+---------+---------+");
        Console.WriteLine("|    N    |   N^2   | Direct  | Recurs. |");
        Console.WriteLine("+---------+---------+---------+---------+");

        // Try it with vectors with the given sizes:
        for(int i=0; i<SIZES.Length; i++) {

            // Computes the average execution time:
            int n = SIZES[i];
            double dtime = Test.TimeIt(FFT.DirectFT, n, REPEAT);
            double rtime = Test.TimeIt(FFT.RecursiveNFFT, n, REPEAT);

            // Print the results:
            string results = String.Format("| {0,7} | {1,7} | {2,7:F4} | {3,7:F4} |",
                n, n*n, dtime, rtime);
            Console.WriteLine(results);
        }

        Console.WriteLine("+---------+---------+---------+---------+");
    }

}
