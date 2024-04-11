/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 **************************************************************************************************
 * This program was compiled and tested in Linux using Mono command line tools. It might need some
 * adaptation to be compiled with Visual Studio (eg. creating a project), but it can be compiles as
 * is with Windows command line tools:
 *
 * $ mcs MainFFT.cs Complex.cs FFT.cs Test.cs
 *
 * This will generate a file called 'MainFFT.exe', that can be run with the command:
 *
 * $ ./MainFFT.exe
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;


class MainFFT
{
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
            double dtime = Test.TimeIt(FFT.DirectFT, n, REPEAT);
            double rtime = Test.TimeIt(FFT.RecursiveFFT, n, REPEAT);
            double itime = Test.TimeIt(FFT.IterativeFFT, n, REPEAT);

            // Print the results:
            string results = String.Format("| {0,7} | {1,7} | {2,7} | {3,7:F4} | {4,7:F4} | {5,7:F4} |",
                n, n*n, n*r, dtime, rtime, itime);
            Console.WriteLine(results);
        }

        Console.WriteLine("+---------+---------+---------+---------+---------+---------+");
    }

}
