/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This library can be compiled with the command:
 *
 * $ javac anyfft.java
 *
 * It can be run by issuing the command:
 *
 * $ java anyfft
 *
 **************************************************************************************************/

/**************************************************************************************************
 * Include necessay libraries:
 **************************************************************************************************/
import java.util.Date;                         // Timing;
import java.util.function.Function;            // Passing functions as parameters;


/**************************************************************************************************
 * Mini-library to deal with complex numbers. While it may be useful to use modules in bigger
 * projects, this one is very small and thorougly tested, and it is simples to just add the
 * code to the main program.
 **************************************************************************************************/
class Complex {

    public float r, i;

    // Constructor:
    public Complex(float re, float im)
    {
        r = re;
        i = im;
    }

    // Constructor:
    public Complex()
    {
        this(0, 0);
    }

    // Add the argument to this, giving the result as a new complex number:
    public Complex add(Complex c)
    {
        return new Complex(r + c.r, i + c.i);
    }

    // Subtract the argument from this, giving the result as a new complex number:
    public Complex sub(Complex c)
    {
        return new Complex(r - c.r, i - c.i);
    }

    // Multiply the argument with this, giving the result as a new complex number:
    public Complex mul(Complex c)
    {
        return new Complex(r*c.r - i*c.i, r*c.i + i*c.r);
    }

    // Divide this by the argument, giving the result as a new complex number:
    public Complex div(float a)
    {
        return new Complex(r/a, i/a);
    }

    // Complex exponential of an angle:
    public static Complex exp(float a)
    {
        return new Complex((float) Math.cos(a), (float) Math.sin(a));
    }

}


/**************************************************************************************************
 * Main class
 **************************************************************************************************/
public class anyfft {

    /**********************************************************************************************
     * Auxiliary Method: complexShow
     *   Pretty printing of an array of complex numbers, used to inspect results.
     *
     * Parameters:
     *   x
     *     A vector of complex numbers, according to the definition above;
     **********************************************************************************************/
    public static void complexShow(Complex x[])
    {
        for(int i=0; i<x.length; i++)
            System.out.printf("( %7.4f, %7.4f )\n", x[i].r, x[i].i);
    }


    /**********************************************************************************************
     * Auxiliary Method: timeIt
     *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
     *
     * Parameters:
     *  f
     *    Function to be called, with the given prototype. The first complex vector is the input
     *    vector, the second complex vector is the result of the computation;
     *  size
     *    Number of elements in the vector on which the transform will be applied;
     *  repeat
     *    Number of times the function will be called.
     *
     * Returns:
     *   The average execution time for that function with a vector of the given size.
     **********************************************************************************************/
    public static float timeIt(Function <Complex[], Complex[]> f, int size, int repeat)
    {
        Complex[] x = new Complex[size];                // Initialize the vector;
        for(int j=0; j<size; j++)
            x[j] = new Complex(j, 0);

        long t0 = System.nanoTime();                    // Start a timer;
        for(int j=0; j<repeat; j++)                     // Repeated calls;
            f.apply(x);
        float ttime = (System.nanoTime() - t0) / (float) (1e9*repeat);
        return ttime;
    }


    /**********************************************************************************************
     * Method: directFT
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
    public static Complex[] directFT(Complex x[])
    {
        int N = x.length;
        Complex[] X = new Complex[N];

        // Initializes twiddle factors;
        Complex W = Complex.exp((float) (-2*Math.PI) / (float) N);
        Complex Wk = new Complex(1, 0);

        for(int k=0; k<N; k++) {
            X[k] = new Complex();                      // Accumulates the results;
            Complex Wkn = new Complex(1, 0);
            for(int n=0; n<N; n++) {
                X[k] = X[k].add(Wkn.mul(x[n]));
                Wkn = Wkn.mul(Wk);                     // Update twiddle factor;
            }
            Wk = Wk.mul(W);
        }
        return X;
    }


    /**********************************************************************************************
     * Method: factor
     *   Smallest prime factor of a given number. If the argument is prime itself, then it is the
     *   return value.
     *
     * Parameters:
     *   n
     *     Number to be inspected.
     *
     * Returns:
     *   The smallest prime factor, or the number itself if it is already a prime.
     **********************************************************************************************/
    public static int factor(int n)
    {
        int rn = n/2;                              // Search up to half the number;
        for(int i=2; i<=rn; i++)
            if (n%i == 0) return i;                // If remainder is zero, a factor is found;
        return n;
    }


    /**********************************************************************************************
     * Method: recursiveFFT
     *   Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
     *   complexity than the direct FT, though the exact value is difficult to compute.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. It must be a composite number, or else the
     *     computation will be defered to the direct FT, and there will be no efficiency gain.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] recursiveFFT(Complex x[])
    {
        int N = x.length;

        int N1 = factor(N);                            // Smallest prime factor of length;
        if (N1==N)                                     // If the length is prime itself,
            return directFT(x);                        //   transform is given by the direct form;
        else {
            int N2 = N / N1;                           // Decompose in two factors, N1 being prime;

            Complex[] xj = new Complex[N2];            // Allocate memory for subsequences;
            Complex[] X = new Complex[N];
            for(int n=0; n<N; n++)
                X[n] = new Complex();

            // Twiddle factors:
            Complex W = Complex.exp((float) (-2*Math.PI) / (float) N);
            Complex Wj = new Complex(1, 0);
            for(int j=0; j<N1; j++) {                  // Compute every subsequence of size N2;
                for(int n=0; n<N2; n++)
                    xj[n] = x[n*N1+j];                 // Create the subsequence;
                Complex[] Xj = recursiveFFT(xj);       // Compute the DFT of the subsequence;
                Complex Wkj = new Complex(1, 0);
                for(int k=0; k<N; k++) {
                    X[k] = X[k].add(Wkj.mul(Xj[k%N2]));        // Recombine results;
                    Wkj = Wkj.mul(Wj);                         // Update twiddle factors;
                }
                Wj = Wj.mul(W);
            }
            return X;
        }
    }


    /**********************************************************************************************
     * Main Method.
     **********************************************************************************************/
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
            double dtime = timeIt(x -> directFT(x), n, REPEAT);
            double rtime = timeIt(x -> recursiveFFT(x), n, REPEAT);

            // Print the results:
            System.out.printf("| %7d | %7d | %7.4f | %7.4f |\n", n, n*n, dtime, rtime);

        }
        System.out.print("+---------+---------+---------+---------+\n");
     }

}
