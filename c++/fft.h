/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This is the implementation file for the FFT.
 **************************************************************************************************/

#ifndef __FFT__
#define __FFT__


// Include necessary libraries:
#include <vector>                              // Deals with arrays;
#include <cmath>                               // Math Functions;
#include "my_complex.h"


using namespace std;


/**************************************************************************************************
 Auxiliary functions:
 **************************************************************************************************/
/**
 * Bit-reversed version of an integer number. This function is private to this module only.
 *
 * @param k The number to be bit-reversed;
 * @param r The number of bits to take into consideration when reversing.
 * @return The number k, bit-reversed according to integers with r bits.
 */
static int bit_reverse(int k, int r)
{
    int l = 0;                                 // Accumulate the results;
    for(int i=0; i<r; i++) {                   // Loop on every bit;
        l = (l << 1) + (k & 1);                // Test less signficant bit and add;
        k >>= 1;                               // Test next bit;
    }
    return l;
}


/**
 * Smallest prime factor of a given number. If the argument is prime itself, then it is the return
 * value. This function is private to this module only.
 *
 * @param n Number to be inspected.
 * @return The smallest prime factor, or the number itself if it is already a prime.
 */
static int factor(int n)
{
    int rn = (int) ceil(sqrt(n));              // Search up to the square root of the number;
    for(int i=2; i<=rn; i++)
        if (n%i == 0) return i;                // If remainder is zero, a factor is found;
    return n;
}


/**************************************************************************************************
 FFT implementations:
 **************************************************************************************************/
/**
 * Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
 *
 * @param x The vector of which the DFT will be computed. Given the nature of the implementation,
 *   there is no restriction on the size of the vector, although it will almost always be called
 *   with a power of two size to give a fair comparison;
 * @return A complex-number vector of the same size, with the coefficients of the DFT.
 */
template<typename T>
vector<Complex<T>> direct_ft(vector<Complex<T>> x)
{
    int N = x.size();
    vector<Complex<T>> X(N);
    Complex<T> W = cexpn<T>(-2*M_PI/N);        // Initialize twiddle factors;
    Complex<T> Wk = Complex<T>(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = Complex<T>();                   // Accumulate the results;
        Complex<T> Wkn = Complex<T>(1, 0);     // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            X[k] = X[k] + Wkn*x[n];
            Wkn = Wkn * Wk;                    // Update twiddle factor;
        }
        Wk = Wk * W;
    }
    return X;
}


/**
 * Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 * complexity.
 *
 * @param x The vector of which the FFT will be computed. This should always be called with a vector
 *   of a power of two length, or it will fail. No checks on this are made.
 * @param[out] X The vector that will receive the results of the computation. It needs to be
 *   allocated prior to the function call;
 * @param N The number of elements in the vector.
 */
template<typename T>
vector<Complex<T>> recursive_fft(vector<Complex<T>> x)
{
    int N = x.size();
    if(N == 1) {                               // A length-1 vector is its own FT;
        vector<Complex<T>> X = x;
        return X;
    } else {
        int N2 = N >> 1;
        vector<Complex<T>> xe(N2);             // Allocate memory for computation;
        vector<Complex<T>> xo(N2);
        vector<Complex<T>> X(N);

        for(int k=0; k<N2; k++) {              // Split even and odd samples;
            xe[k] = x[k<<1];
            xo[k] = x[(k<<1)+1];
        }
        vector<Complex<T>> Xe = recursive_fft(xe);     // Transform of even samples;
        vector<Complex<T>> Xo = recursive_fft(xo);     // Transform of odd samples;

        Complex<T> W = cexpn<T>(-2*M_PI/N);            // Twiddle factors;
        Complex<T> Wk = Complex<T>(1, 0);
        for(int k=0; k<N2; k++) {
            Complex<T> WkXo = Wk * Xo[k];              // Recombine results;
            X[k] = Xe[k] + WkXo;
            X[k+N2] = Xe[k] - WkXo;
            Wk = Wk * W;                               // Update twiddle factors;
        }
        return X;
    }
}


/**
 * Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
 * O(N log_2(N)) complexity, and since there are less function calls, it will probably be marginally
 * faster than the recursive versions.
 *
 * @param x The vector of which the FFT will be computed. This should always be called with a vector
 *   of a power of two length, or it will fail. No checks on this are made.
 * @param[out] X The vector that will receive the results of the computation. It needs to be
 *   allocated prior to the function call;
 * @param N The number of elements in the vector.
 */
template<typename T>
vector<Complex<T>> iterative_fft(vector<Complex<T>> x)
{
    int N = x.size();
    int r = (int) floor(log2(N));              // Number of bits;
    vector<Complex<T>> X(N);
    for(int k=0; k<N; k++) {
        int l = bit_reverse(k, r);             // Reorder the vector according to the
        X[l] = x[k];                           //   bit-reversed order;
    }

    int step = 1;                              // Auxiliary for computation of twiddle factors;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            Complex<T> W = cexpn<T>(-M_PI/step);       // Twiddle factors;
            Complex<T> Wkn = Complex<T>(1, 0);
            for(int n=0; n<step; n++) {
                int p = l + n;
                int q = p + step;
                X[q] = X[p] - Wkn * X[q];              // Recombine results;
                X[p] = X[p]*2 - X[q];
                Wkn = Wkn * W;                         // Update twiddle factors;
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
 * @param x The vector of which the FFT will be computed. Its length must be a composite number, or
 *   else the computation will be defered to the direct FT, and there will be no efficiency gain.
 * @param[out] X The vector that will receive the results of the computation. It needs to be
 *   allocated prior to the function call;
 * @param N The number of elements in the vector.
 */
template<typename T>
vector<Complex<T>> recursive_nfft(vector<Complex<T>> x)
{
    int N = x.size();
    int N1 = factor(N);                                // Smallest prime factor of length;
    if(N1 == N) {                                      // If the length is prime itself,
        vector<Complex<T>> X = direct_ft<T>(x);        //   transform is given by the direct form;
        return X;
    } else {
        int N2 = N / N1;                               // Decompose in two factors, N1 being prime;
        vector<Complex<T>> xj(N2);                     // Allocate memory for subsequences;
        vector<Complex<T>> X(N);

        Complex<T> W = cexpn<T>(-2*M_PI/N);            // Twiddle factor;
        Complex<T> Wj = Complex<T>(1, 0);
        for(int j=0; j<N1; j++) {                      // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++)
                xj[n] = x[n*N1+j];                     // Create the subsequence;
            vector<Complex<T>> Xj = recursive_nfft<T>(xj);     // Compute the DFT of the subsequence;
            Complex<T> Wkj = Complex<T>(1, 0);
            for(int k=0; k<N; k++) {
                X[k] = X[k] + Xj[k%N2] * Wkj;          // Recombine results;
                Wkj = Wkj * Wj;                        // Update twiddle factors;
            }
            Wj = Wj * W;
        }
        return X;
    }
}


#endif
