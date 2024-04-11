/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This version implements Cooley-Tukey algorithm.
 **************************************************************************************************/

package fft;


public class FFT {

    /**
     * Bit-reversed version of an integer number.
     *
     * @param k The number to be bit-reversed;
     * @param r The number of bits to take into consideration when reversing.
     * @return The number k, bit-reversed according to integers with r bits.
     */
    private static int bitReverse(int k, int r)
    {
        int l = 0;                                     // Accumulate the results;
        for(int i=0; i<r; i++) {                       // Loop on every bit;
            l = (l << 1) + (k & 1);                    // Test less signficant bit and add;
            k = (k >> 1);                              // Test next bit;
        }
        return l;
    }


    /**
     * Smallest prime factor of a given number. If the argument is prime itself, then it is the
     * return value.
     *
     * @param n Number to be inspected.
     * @return The smallest prime factor, or the number itself if it is already a prime.
     */
    public static int factor(int n)
    {
        int rn = (int) Math.ceil(Math.sqrt(n));    // Search up to the square root of the number;
        for(int i=2; i<=rn; i++)
            if (n%i == 0) return i;                // If remainder is zero, a factor is found;
        return n;
    }


    /**
     * Discrete Fourier Transform directly from the definition. This algorithm has O(N^2)
     * complexity.
     *
     * @param x The vector of which the DFT will be computed. Given the nature of the
     *          implementation, there is no restriction on the size of the vector, although it will
     *          almost always be called with a power of two size to give a fair comparison;
     * @return A complex-number vector of the same size, with the coefficients of the DFT.
     */
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


    /**
     * Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
     * complexity.
     *
     * @param x The vector of which the FFT will be computed. This should always be called with a
     *          vector of a power of two length, or it will fail. No checks on this are made.
     * @return A complex-number vector of the same size, with the coefficients of the DFT.
     */
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


    /**
     * Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
     * O(N log_2(N)) complexity, and since there are less function calls, it will probably be
     * marginally faster than the recursive versions.
     *
     * @param x The vector of which the FFT will be computed. This should always be called with a
     *          vector of a power of two length, or it will fail. No checks on this are made.
     * @return A complex-number vector of the same size, with the coefficients of the DFT.
     */
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


    /**
     * Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
     * complexity than the direct FT, though the exact value is difficult to compute.
     *
     * @param x The vector of which the FFT will be computed. Its length must be a composite number,
     *          or else the computation will be defered to the direct FT, and there will be no
     *          efficiency gain.
     * @return A complex-number vector of the same size, with the coefficients of the DFT.
     */
    public static Complex[] recursiveNFFT(Complex x[])
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
                Complex[] Xj = recursiveNFFT(xj);      // Compute the DFT of the subsequence;
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

}
