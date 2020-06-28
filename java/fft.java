/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This library can be compiled with the command:
 *
 * $ javac fft.java
 *
 * It can be run by issuing the command:
 *
 * $ java fft
 *
 **************************************************************************************************/

/**************************************************************************************************
 * Include necessary libraries:
 **************************************************************************************************/
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
public class fft {

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
        int N = x.length;                              // Length of the vector;
        Complex[] X = new Complex[N];                  // Accumulate the results;

        // Initializes twiddle factors:
        Complex W = Complex.exp((float) (-2*Math.PI) / (float) N);
        Complex Wk = new Complex(1, 0);

        for(int k=0; k<N; k++) {                       // Compute the kth coefficient;
            X[k] = new Complex();                      // Accumulate the results;
            Complex Wkn = new Complex(1, 0);
            for(int n=0; n<N; n++) {                   //   Operate the summation;
                X[k] = X[k].add(Wkn.mul(x[n]));        //     Compute every term;
                Wkn = Wkn.mul(Wk);                     // Update twiddle factor;
            }
            Wk = Wk.mul(W);
        }
        return X;
    }


    /**********************************************************************************************
     * Method: recursiveFFT
     *   Fast Fourier Transform using a recursive decimation in time algorithm. This has
     *   O(N log_2(N)) complexity.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. This should always be called with a vector
     *     of a power of two length, or it will fail. No checks on this are made.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] recursiveFFT(Complex x[])
    {
        int N = x.length;

        if(N==1)                                       // A length-1 vector is its own FT;
            return x;
        else {
            int N2 = N >> 1;

            Complex[] xe = new Complex[N2];            // Allocate memory for computation;
            Complex[] xo = new Complex[N2];
            Complex[] X = new Complex[N];

            for(int k=0; k<N2; k++) {                  // Split even and odd samples;
                xe[k] = x[k<<1];
                xo[k] = x[(k<<1)+1];
            }
            Complex[] Xe = recursiveFFT(xe);           // Transform of even samples;
            Complex[] Xo = recursiveFFT(xo);           // Transform of odd samples;

            // Twiddle factors:
            Complex W = Complex.exp((float) (-2*Math.PI) / (float) N);
            Complex Wk = new Complex(1, 0);
            for(int k=0; k<N2; k++) {
                Complex w = Wk.mul(Xo[k]);             // Recombine results;
                X[k] = Xe[k].add(w);
                X[k+N2] = Xe[k].sub(w);
                Wk = Wk.mul(W);                        // Update twiddle factors;
            }
            return X;
        }
    }


    /**********************************************************************************************
     * Method: bitReverse
     *   Bit-reversed version of an integer number.
     *
     * Parameters:
     *   k
     *     The number to be bit-reversed;
     *   r
     *     The number of bits to take into consideration when reversing.
     *
     * Returns:
     *   The number k, bit-reversed according to integers with r bits.
     **********************************************************************************************/
    private static int bitReverse(int k, int r)
    {
        int l = 0;                                     // Accumulate the results;
        for(int i=0; i<r; i++) {                       // Loop on every bit;
            l = (l << 1) + (k & 1);                    // Test less signficant bit and add;
            k = (k >> 1);                              // Test next bit;
        }
        return l;
    }


    /**********************************************************************************************
     * Function: iterativeFFT
     *   Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
     *   O(N log_2(N)) complexity, and since there are less function calls, it will probably be
     *   marginally faster than the recursive versions.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. This should always be called with a vector
     *     of a power of two length, or it will fail. No checks on this are made.
     *
     * Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    public static Complex[] iterativeFFT(Complex x[])
    {
        int N = x.length;
        Complex[] X = new Complex[N];

        int r = (int) Math.round(Math.log(N)/Math.log(2));     // Number of bits;
        for(int k=0; k<N; k++) {
            int l = bitReverse(k, r);                          // Reorder the vector according to
            X[l] = x[k];                                       //   the bit-reversed order;
        }

        int step = 1;                                          // Computation of twiddle factors;
        for(int k=0; k<r; k++) {
            for(int l=0; l<N; l+=2*step) {
                // Twiddle factors:
                Complex W = Complex.exp((float) (-Math.PI) / (float) step);
                Complex Wkn = new Complex(1, 0);
                for(int n=0; n<step; n++) {
                    int p = l + n;
                    int q = p + step;
                    X[q] = X[p].sub(Wkn.mul(X[q]));            // Recombine results;
                    X[p] = (new Complex(2,0)).mul(X[p]).sub(X[q]);
                    Wkn = Wkn.mul(W);                          // Update twiddle factors;
                }
            }
            step <<= 1;
        }
        return X;
    }

    /**********************************************************************************************
     * Main Method.
     **********************************************************************************************/
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
            double dtime = timeIt(x -> directFT(x), n, REPEAT);
            double rtime = timeIt(x -> recursiveFFT(x), n, REPEAT);
            double itime = timeIt(x -> iterativeFFT(x), n, REPEAT);

            // Print the results:
            System.out.printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                              n, n*n, n*r, dtime, rtime, itime);

        }
        System.out.print("+---------+---------+---------+---------+---------+---------+\n");
    }

}
