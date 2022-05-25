/**************************************************************************************************
 * Fast Fourier Transform -- C with Native Complex Version
 * This is the implementation file for the FFT using native C complex numbers.
 **************************************************************************************************/

// Include necessary libraries:
#include "fft.h"
#include "fft_float.h"


void float_direct_ft(float complex x[], float complex X[], int N)
{
    float complex W, Wk, Wkn;                  // Twiddle factors;

    W = cexpf(-2i*M_PI/N);                     // Initialize twiddle factors;
    Wk = 1;
    for(int k=0; k<N; k++) {
        X[k] = 0;                              // Accumulate the results;
        Wkn = 1;                               // Initialize twiddle factors;
        for(int n=0; n<N; n++) {
            X[k] = X[k] + Wkn * x[n];
            Wkn = Wkn * Wk;                    // Update twiddle factors;
        }
        Wk = Wk * W;
    }
}


void float_recursive_fft(float complex x[], float complex X[], int N)
{
    float complex *xe, *xo, *Xe, *Xo;          // Vectors with intermediate results;
    float complex W, Wk, WkXo;                 // Twiddle factors;
    int N2;

    if(N==1)                                   // A length-1 vector is its own FT;
        X[0] = x[0];
    else {
        N2 = N >> 1;

        xe = (float complex *) malloc(N2 * sizeof(float complex));     // Allocate memory;
        xo = (float complex *) malloc(N2 * sizeof(float complex));
        Xe = (float complex *) malloc(N2 * sizeof(float complex));
        Xo = (float complex *) malloc(N2 * sizeof(float complex));

        for(int k=0; k<N2; k++) {                      // Split even and odd samples;
            xe[k] = x[k<<1];
            xo[k] = x[(k<<1)+1];
        }
        float_recursive_fft(xe, Xe, N2);               // Transform of even samples;
        float_recursive_fft(xo, Xo, N2);               // Transform of odd samples;

        W = cexpf(-2i*M_PI/N);                         // Twiddle factors;
        Wk = 1;
        for(int k=0; k<N2; k++) {
            WkXo = Wk * Xo[k];
            X[k] = Xe[k] + WkXo;
            X[k+N2] = Xe[k] - WkXo;
            Wk = Wk * W;                               // Update twiddle factors;
        }

        free(Xo);                                      // Releasing memory of intermediate vectors;
        free(Xe);
        free(xo);
        free(xe);
    }
}


void float_iterative_fft(float complex x[], float complex X[], int N)
{
    float complex W, Wkn;                      // Twiddle factors;
    int r, l, p, q, step;

    r = (int) floor(log2(N));                  // Number of bits;
    for(int k=0; k<N; k++) {
        l = bit_reverse(k, r);                 // Reorder the vector according to the
        X[l] = x[k];                           //   bit-reversed order;
    }

    step = 1;                                  // Auxiliary for computation of twiddle factors;
    for(int k=0; k<r; k++) {
        for(int l=0; l<N; l+=2*step) {
            W = cexpf(-1i*M_PI/step);          // Twiddle factors;
            Wkn = 1;
            for(int n=0; n<step; n++) {
                p = l + n;
                q = p + step;
                X[q] = X[p] - Wkn * X[q];
                X[p] = 2 * X[p] - X[q];
                Wkn = Wkn * W;                 // Update twiddle factors;
            }
        }
        step <<= 1;
    }
}


void float_recursive_nfft(float complex x[], float complex X[], int N)
{
    float complex *xj, *Xj;                    // Subsequences of the transform;
    float complex W, Wj, Wkj;
    int N1, N2;

    N1 = factor(N);                            // Smallest prime factor of the length;
    if (N1==N)                                 // If the length is prime itself,
        float_direct_ft(x, X, N);              //   the transform is given by the direct form;
    else {
        N2 = N / N1;                           // Decompose in two factors, N1 being prime;

        xj = malloc(N2 * sizeof(float complex));       // Allocate memory for subsequences
        Xj = malloc(N2 * sizeof(float complex));       //    and their transforms;

        W = cexpf(-2i*M_PI/N);                 // Twiddle factor;
        Wj = 1;
        for(int j=0; j<N1; j++) {              // Compute every subsequence of size N2;
            for(int n=0; n<N2; n++) {
                xj[n] = x[n*N1+j];             // Create the subsequence;
                Xj[n] = 0;                     // Initialize the transform;
            }
            float_recursive_nfft(xj, Xj, N2);  // Compute the DFT of the subsequence;
            Wkj = 1;
            for(int k=0; k<N; k++) {
                X[k] = X[k] + Wkj * Xj[k%N2];  // Recombine results;
                Wkj = Wkj * Wj;                // Update twiddle factors;
            }
            Wj = Wj * W;
        }

        free(xj);                              // Clean-up used memory;
        free(Xj);
    }
}
