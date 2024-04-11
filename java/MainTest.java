/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This program performs a simple test of the implementations.
 **************************************************************************************************
 * This program can be compiled with the command:
 *
 * $ javac -cp . MainTest.java
 *
 * It can be run by issuing the command:
 *
 * $ java MainTest
 **************************************************************************************************/

// Include necessary libraries:
import java.util.function.Function;            // Passing functions as parameters;
import fft.*;                                  // FFT Package;


class MainTest {

    public static void main(String args[])
    {
        // Tests for the implementations of the DFT of power-of-2 length:
        System.out.println("Direct FT - ");
        Test.testIt(x -> FFT.directFT(x), 8);
        System.out.println("Recursive FFT - ");
        Test.testIt(x -> FFT.recursiveFFT(x), 8);
        System.out.println("Iterative FFT - ");
        Test.testIt(x -> FFT.iterativeFFT(x), 8);
        System.out.println("Direct FT - ");
        Test.testIt(x -> FFT.directFT(x), 16);
        System.out.println("Recursive FFT - ");
        Test.testIt(x -> FFT.recursiveFFT(x), 16);
        System.out.println("Iterative FFT - ");
        Test.testIt(x -> FFT.iterativeFFT(x), 16);

        // Tests for the implementations of the DFT of composite length:
        System.out.println("Direct FT - ");
        Test.testIt(x -> FFT.directFT(x), 12);
        System.out.println("Recursive FFT - ");
        Test.testIt(x -> FFT.recursiveNFFT(x), 12);
    }

}
