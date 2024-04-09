/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 **************************************************************************************************
 * This program can be compiled with the command:
 *
 * $ javac -cp . MainAnyFFT.java
 *
 * It can be run by issuing the command:
 *
 * $ java MainAnyFFT
 **************************************************************************************************/

// Include necessary libraries:
import java.util.function.Function;            // Passing functions as parameters;
import fft.*;                                  // FFT Package;


public class MainAnyFFT {

    public static void main(String args[])
    {
        int[] SIZES = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
        int REPEAT = 500;                      // Number of executions to compute average time;

        // Start by printing the table with time comparisons:
        System.out.print("+---------+---------+---------+---------+\n");
        System.out.print("|    N    |   N^2   | Direct  | Recurs. |\n");
        System.out.print("+---------+---------+---------+---------+\n");

        // Try it with vectors with the given sizes:
        for(int i=0; i<SIZES.length; i++) {

            // Compute the average execution time:
            int n = SIZES[i];
            double dtime = Test.timeIt(x -> FFT.directFT(x), n, REPEAT);
            double rtime = Test.timeIt(x -> FFT.recursiveNFFT(x), n, REPEAT);

            // Print the results:
            System.out.printf("| %7d | %7d | %7.4f | %7.4f |\n", n, n*n, dtime, rtime);

        }
        System.out.print("+---------+---------+---------+---------+\n");
     }

}
