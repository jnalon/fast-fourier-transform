/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This version implements Cooley-Tukey algorithm.
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;


class FFT
{
    // Prototype to the DFT function, to be used at the timing function:
    public delegate Complex[] DFT(Complex[] x);


    /// <summary>
    ///  Bit-reversed version of an integer number.
    /// <param name="k">The number to be bit-reversed.</param>
    /// <param name="r">The number of bits to take into consideration when reversing.</param>
    /// <returns>The number k, bit-reversed according to integers with r bits.</returns>
    /// </summary>
    private static int BitReverse(int k, int r)
    {
        int l = 0;                                     // Accumulate the results;
        for(int i=0; i<r; i++) {                       // Loop on every bit;
            l = (l << 1) + (k & 1);                    // Test less signficant bit and add;
            k = (k >> 1);                              // Test next bit;
        }
        return l;
    }


    /// <summary>
    /// Smallest prime factor of a given number. If the argument is prime itself, then it is the
    /// return value.
    /// <param name="n">Number to be inspected.</param>
    /// <returns>
    ///    The smallest prime factor, or the number itself if it is already a prime.
    /// </returns>
    /// </summary>
    public static int Factor(int n)
    {
    int rn = (int) Math.Ceiling(Math.Sqrt(n));     // Search up to the square root of the number;
        for(int i=2; i<=rn; i++)
            if (n%i == 0) return i;                // If remainder is zero, a factor is found;
        return n;
    }


    /// <summary>
    /// Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    /// complexity.
    /// <param name="x">
    ///    The vector of which the DFT will be computed. Given the nature of the implementation,
    ///    there is no restriction on the size of the vector, although it will almost always be
    ///    called with a power of two size to give a fair comparison.
    /// </param>
    /// <returns>
    ///     A complex-number vector of the same size, with the coefficients of the DFT.
    /// </returns>
    /// </summary>
    public static Complex[] DirectFT(Complex[] x)
    {
        int N = x.Length;                              // Length of the vector;
        Complex[] X = new Complex[N];                  // Accumulate the results;
        Complex W = Complex.exp(-2*Math.PI/N);         // Initialize twiddle factors;
        Complex Wk = new Complex(1, 0);

        for(int k=0; k<N; k++) {                       // Compute the kth coefficient;
            X[k] = new Complex();                      // Accumulate the results;
            Complex Wkn = new Complex(1, 0);
            for(int n=0; n<N; n++) {                   //   Operate the summation;
                X[k] = X[k] + Wkn*x[n];                //     Compute every term;
                Wkn = Wkn * Wk;                        // Update twiddle factor;
            }
            Wk = Wk * W;
        }
        return X;
    }


    /// <summary>
    /// Fast Fourier Transform using a recursive decimation in time algorithm. This has
    /// O(N log_2(N)) complexity.
    /// <param name="x">
    ///    The vector of which the FFT will be computed. This should always be called with a vector
    ///    of a power of two length, or it will fail. No checks on this are made.
    /// </param>
    /// <returns>
    ///  A complex-number vector of the same size, with the coefficients of the DFT.
    /// </returns>
    /// </summary>
    public static Complex[] RecursiveFFT(Complex[] x)
    {
        int N = x.Length;

        if (N==1)                                      // A length-1 vector is its own FT;
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
            Complex[] Xe = RecursiveFFT(xe);           // Transform of even samples;
            Complex[] Xo = RecursiveFFT(xo);           // Transform of odd samples;

            Complex W = Complex.exp(-2*Math.PI/N);     // Twiddle factors;
            Complex Wk = new Complex(1, 0);
            for(int k=0; k<N2; k++) {
                Complex w = Wk * Xo[k];                // Recombine results;
                X[k] = Xe[k] + w;
                X[k+N2] = Xe[k] - w;
                Wk = Wk * W;                           // Update twiddle factors;
            }
            return X;
        }
    }


    /// <summary>
    /// Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
    /// O(N log_2(N)) complexity, and since there are less function calls, it will probably be
    /// marginally faster than the recursive versions.
    /// <param name="x">
    ///    The vector of which the FFT will be computed. This should always be called with a vector
    ///    of a power of two length, or it will fail. No checks on this are made.
    /// </param>
    /// <returns>
    ///    A complex-number vector of the same size, with the coefficients of the DFT.
    /// </returns>
    /// </summary>
    public static Complex[] IterativeFFT(Complex[] x)
    {
        int N = x.Length;
        Complex[] X = new Complex[N];

        int r = (int) Math.Round(Math.Log(N, 2));              // Number of bits;
        for(int k=0; k<N; k++) {
            int l = BitReverse(k, r);                          // Reorder the vector according to
            X[l] = x[k];                                       //   the bit-reversed order;
        }

        int step = 1;                                          // Computation of twiddle factors;
        for(int k=0; k<r; k++) {
            for(int l=0; l<N; l+=2*step) {
                Complex W = Complex.exp(-Math.PI/step);        // Twiddle factors;
                Complex Wkn = new Complex(1, 0);
                for(int n=0; n<step; n++) {
                    int p = l + n;
                    int q = p + step;
                    X[q] = X[p] - Wkn * X[q];                  // Recombine results;
                    X[p] = X[p]*2 - X[q];
                    Wkn = Wkn * W;                             // Update twiddle factors;
                }
            }
            step <<= 1;
        }

        return X;
    }


    /// <summary>
    /// Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
    /// complexity than the direct FT, though the exact value is difficult to compute.
    /// <param name="x">
    ///    The vector of which the FFT will be computed. Its length must be a composite number, or
    ///    else the computation will be defered to the direct FT, and there will be no efficiency
    ///    gain.
    /// </param>
    /// <returns>
    ///    A complex-number vector of the same size, with the coefficients of the DFT.
    /// </returns>
    /// </summary>
    public static Complex[] RecursiveNFFT(Complex[] x)
    {
        int N = x.Length;
        Complex[] X = new Complex[N];
        for(int n=0; n<N; n++)
            X[n] = new Complex();

        int N1 = Factor(N);                            // Smallest prime factor of length;
        if (N1==N)                                     // If the length is prime itself,
            return DirectFT(x);                        //   transform is given by the direct form;
        else {
            int N2 = N / N1;                           // Decompose in two factors, N1 being prime;

            Complex[] xj = new Complex[N2];            // Allocate memory for subsequences;

            Complex W = Complex.exp(-2*Math.PI/N);     // Twiddle factor;
            Complex Wj = new Complex(1, 0);
            for(int j=0; j<N1; j++) {                  // Compute every subsequence of size N2;
                for(int n=0; n<N2; n++)
                    xj[n] = x[n*N1+j];                 // Create the subsequence;
                Complex[] Xj = RecursiveNFFT(xj);      // Compute the DFT of the subsequence;
                Complex Wkj = new Complex(1, 0);
                for(int k=0; k<N; k++) {
                    X[k] = X[k] + Xj[k%N2] * Wkj;      // Recombine results;
                    Wkj = Wkj * Wj;                    // Update twiddle factors;
                }
                Wj = Wj * W;
            }
            return X;
        }
    }


}
