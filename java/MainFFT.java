/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 **************************************************************************************************
 * This program can be compiled with the command:
 *
 * $ javac -cp . MainFFT.java
 *
 * It can be run by issuing the command:
 *
 * $ java MainFFT
 **************************************************************************************************/

// Include necessary libraries:
import java.util.function.Function;            // Passing functions as parameters;
import fft.*;                                  // FFT Package;


public class MainFFT {

    public static void main(String args[])
    {
        int REPEAT = 500;                      // Number of executions to compute average time;

        // Start by printing the table with time comparisons:
        System.out.print("+---------+---------+---------+---------+---------+---------+\n");
        System.out.print("|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |\n");
        System.out.print("+---------+---------+---------+---------+---------+---------+\n");

        // Try it with vectors with size ranging from 32 to 1024 samples:
        for(int r=5; r<11; r++) {

            // Compute the average execution time:
            int n = (int) Math.pow(2, r);
            double dtime = Test.timeIt(x -> FFT.directFT(x), n, REPEAT);
            double rtime = Test.timeIt(x -> FFT.recursiveFFT(x), n, REPEAT);
            double itime = Test.timeIt(x -> FFT.iterativeFFT(x), n, REPEAT);

            // Print the results:
            System.out.printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                              n, n*n, n*r, dtime, rtime, itime);

        }
        System.out.print("+---------+---------+---------+---------+---------+---------+\n");
    }

}
