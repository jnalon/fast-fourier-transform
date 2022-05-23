/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the implementation file for the FFT.
 **************************************************************************************************/

// Include necessary libraries:
#include "my_complex.h"
#include "fft.h"


void direct_ft(Complex x[], Complex X[], int N)
{
    Complex W, Wk, Wkn, w;                     // Twiddle factors;

    W = cexpn(-2*M_PI/N);                      // Initialize twiddle factors;
    Wk = cmplx(1, 0);
    for(int k=0; k<N; k++) {
        X[k] = cmplx(0, 0);                    // Accumulate the results;
        Wkn = cmplx(1, 0);                     // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            w = cmul(x[n], Wkn);
            X[k] = cadd(X[k], w);
            Wkn = cmul(Wkn, Wk);               // Update twiddle factors;
        }
        Wk = cmul(Wk, W);
    }
}


void recursive_fft(Complex x[], Complex X[], int N)
{
    Complex *xe, *xo, *Xe, *Xo;                // Vectors with intermediate results;
    Complex W, Wk, WkXo;                       // Twiddle factors;
    int N2;

    if(N==1)                                   // A length-1 vector is its own FT;
        X[0] = x[0];
    else {
        N2 = N >> 1;

        xe = (Complex *) malloc(N2 * sizeof(Complex)); // Allocate memory for computation;
        xo = (Complex *) malloc(N2 * sizeof(Complex));
        Xe = (Complex *) malloc(N2 * sizeof(Complex));
        Xo = (Complex *) malloc(N2 * sizeof(Complex));

        for(int k=0; k<N2; k++) {                      // Split even and odd samples;
            xe[k] = x[k<<1];
            xo[k] = x[(k<<1)+1];
            Xe[k] = cmplx(0, 0);
            Xo[k] = cmplx(0, 0);
        }
        recursive_fft(xe, Xe, N2);                     // Transform of even samples;
        recursive_fft(xo, Xo, N2);                     // Transform of odd samples;

        W = cexpn(-2*M_PI/N);                          // Twiddle factors;
        Wk = cmplx(1, 0);
        for(int k=0; k<N2; k++) {
            WkXo = cmul(Xo[k], Wk);
            X[k] = cadd(Xe[k], WkXo);                  // Recombine results;
            X[k+N2] = csub(Xe[k], WkXo);
            Wk = cmul(Wk, W);                          // Update twiddle factors;
        }

        free(Xo);                                      // Releasing memory of intermediate vectors;
        free(Xe);
        free(xo);
        free(xe);
    }
}


int bit_reverse(int k, int r)
{
    int l = 0;                                 // Accumulate the results;
    for(int i=0; i<r; i++) {                   // Loop on every bit;
        l = (l << 1) + (k & 1);                // Test less signficant bit and add;
        k >>= 1;                               // Test next bit;
    }
    return l;
}


void iterative_fft(Complex x[], Complex X[], int N)
{
    Complex W, Wkn, w;                         // Twiddle factors;
    int r, l, p, q, step;

    r = (int) floor(log2(N));                  // Number of bits;
    for(int k=0; k<N; k++) {
        l = bit_reverse(k, r);                 // Reorder the vector according to the
        X[l] = x[k];                           //   bit-reversed order;
    }

    step = 1;                                  // Auxiliary for computation of twiddle factors;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            W = cexpn(-M_PI/step);             // Twiddle factors;
            Wkn = cmplx(1, 0);
            for(int n=0; n<step; n++) {
                p = l + n;
                q = p + step;
                w = cmul(X[q], Wkn);           // Recombine results;
                X[q] = csub(X[p], w);
                w = cmul(X[p], cmplx(2, 0));
                X[p] = csub(w, X[q]);
                Wkn = cmul(Wkn, W);            // Update twiddle factors;
            }
        }
        step <<= 1;
    }
}


int factor(int n)
{
    int rn = (int) ceil(sqrt(n));              // Search up to the square root of the number;
    for(int i=2; i<=rn; i++)
        if (n%i == 0) return i;                // If remainder is zero, a factor is found;
    return n;
}


void recursive_nfft(Complex x[], Complex X[], int N)
{
    Complex *xj, *Xj;                          // Subsequences of the transform;
    Complex W, Wj, Wkj;
    int N1, N2;

    N1 = factor(N);                            // Smallest prime factor of the length;
    if (N1==N)                                 // If the length is prime itself,
        direct_ft(x, X, N);                    //   the transform is given by the direct form;
    else {
        N2 = N / N1;                           // Decompose in two factors, N1 being prime;

        xj = malloc(N2 * sizeof(Complex));     // Allocate memory for subsequences
        Xj = malloc(N2 * sizeof(Complex));     //    and their transforms;

        W = cexpn(-2*M_PI/N);                  // Twiddle factor;
        Wj = cmplx(1, 0);
        for(int j=0; j<N1; j++) {              // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++) {
                xj[n] = x[n*N1+j];             // Create the subsequence;
                Xj[n] = cmplx(0, 0);           // Initialize the transform;
            }
            recursive_nfft(xj, Xj, N2);        // Compute the DFT of the subsequence;
            Wkj = cmplx(1, 0);
            for(int k=0; k<N; k++) {
                X[k] = cadd(X[k], cmul(Xj[k%N2], Wkj));        // Recombine results;
                Wkj = cmul(Wkj, Wj);                           // Update twiddle factors;
            }
            Wj = cmul(Wj, W);
        }
        free(xj);                              // Clean-up used memory;
        free(Xj);
    }
}
