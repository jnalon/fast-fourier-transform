/**************************************************************************************************
 * Fast Fourier Transform -- C with Native Complex Version
 * This is the header file for the implementations of the FFT using C native complex numbers.
 **************************************************************************************************/

#ifndef __FFT_DOUBLE__
#define __FFT_DOUBLE__


// Include necessary libraries:
#include <stdlib.h>                            // Standard Library;
#include <stdio.h>                             // Input and Output;
#include <math.h>                              // Math Functions;
#include <complex.h>                           // Complex numbers;


/**
 * Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
 *
 * @param x The vector of which the DFT will be computed. Given the nature of the implementation,
 *   there is no restriction on the size of the vector, although it will almost always be called
 *   with a power of two size to give a fair comparison;
 * @param[out] X The vector that will receive the results of the computation. It needs to be
 *   allocated prior to the function call;
 * @param N The number of elements in the vector.
 */
void double_direct_ft(double complex x[], double complex X[], int N);


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
void double_recursive_fft(double complex x[], double complex X[], int N);


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
void double_iterative_fft(double complex x[], double complex X[], int N);


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
void double_recursive_nfft(double complex x[], double complex X[], int N);


#endif
