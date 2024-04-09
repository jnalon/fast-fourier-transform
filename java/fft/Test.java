/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This file implements methods for testing and time measuring.
 **************************************************************************************************/

package fft;


// Include necessary libraries:
import java.util.function.Function;            // Passing functions as parameters;


public class Test {

    /**
     * Initializes a complex vector. This method creates and initializes a vector of complex numbers
     * with bogus data for testing purposes only.
     *
     * @param size Number of elements in the vector.
     * @return The initialized vector.
     */
    private static Complex[] initializeVector(int size)
    {
        Complex[] x = new Complex[size];                // Initialize the vector;
        for(int j=0; j<size; j++)
            x[j] = new Complex(j, 0);
        return x;
    }


    /**
     * Pretty printing of an array of complex numbers, used to inspect results.
     *
     * @param f Function to be called, with the given prototype. The first complex vector is the
     *          input vector, the second complex vector is the result of the computation;
     * @param size Number of elements in the vector on which the transform will be applied;
     */
    public static void testIt(Function <Complex[], Complex[]> f, int size)
    {
        Complex[] x = initializeVector(size);
        Complex[] X = f.apply(x);
        System.out.printf("N = %d | Input | Output:\n", size);
        for(int i=0; i<size; i++)
            System.out.printf("  %2d: ( %8.4f, %8.4f ) | ( %8.4f, %8.4f )\n",
                              i, x[i].r, x[i].i, X[i].r, X[i].i);
        System.out.println("------------------------------");
    }


    /**
     * Measure execution time through repeated calls to a (Fast) Fourier Transform function.
     *
     * @param f Function to be called, with the given prototype. The first complex vector is the
     *          input vector, the second complex vector is the result of the computation;
     * @param size Number of elements in the vector on which the transform will be applied;
     * @param repeat Number of times the function will be called.
     * @return The average execution time for that function with a vector of the given size.
     */
    public static float timeIt(Function <Complex[], Complex[]> f, int size, int repeat)
    {
        Complex[] x = initializeVector(size);
        long t0 = System.nanoTime();                    // Start a timer;
        for(int j=0; j<repeat; j++)                     // Repeated calls;
            f.apply(x);
        float ttime = (System.nanoTime() - t0) / (float) (1e9*repeat);
        return ttime;
    }

}
