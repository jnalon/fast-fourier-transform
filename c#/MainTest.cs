/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This program performs a simple test of the implementations.
 **************************************************************************************************
 * This program was compiled and tested in Linux using Mono command line tools. It might need some
 * adaptation to be compiled with Visual Studio (eg. creating a project), but it can be compiles as
 * is with Windows command line tools:
 *
 * $ mcs MainTest.cs FFT.cs Test.cs Complex.cs
 *
 * This will generate a file called 'MainTest.exe', that can be run with the command:
 *
 * $ ./MainFFT.exe
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;
using System.Diagnostics;                      // Time measuring;


class MainTest
{

    public static void Main()
    {
        // Tests for the implementations of the DFT of power-of-2 length:
        Console.WriteLine("Direct FT - ");
        Test.TestIt(FFT.DirectFT, 8);
        Console.WriteLine("Recursive FFT - ");
        Test.TestIt(FFT.RecursiveFFT, 8);
        Console.WriteLine("Iterative FFT - ");
        Test.TestIt(FFT.IterativeFFT, 8);
        Console.WriteLine("Direct FT - ");
        Test.TestIt(FFT.DirectFT, 16);
        Console.WriteLine("Recursive FFT - ");
        Test.TestIt(FFT.RecursiveFFT, 16);
        Console.WriteLine("Iterative FFT - ");
        Test.TestIt(FFT.IterativeFFT, 16);

        // Tests for the implementations of the DFT of composite length:
        Console.WriteLine("Direct FT - ");
        Test.TestIt(FFT.DirectFT, 12);
        Console.WriteLine("Recursive FFT - ");
        Test.TestIt(FFT.RecursiveNFFT, 12);
    }


}
