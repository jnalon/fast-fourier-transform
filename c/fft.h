/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the header file for the implementations of the FFT.
 **************************************************************************************************/

#ifndef __FFT__
#define __FFT__


// Include necessary libraries:
#include <stdlib.h>                            // Standard Library;
#include <math.h>                              // Math Functions;
#include "my_complex.h"                        // Mini-library for complex numbers;


/**************************************************************************************************
 Auxiliary functions:
 **************************************************************************************************/
/**
 * Bit-reversed version of an integer number. This function is not exported to the header, as it is
 * only used internally
 *
 * @param k The number to be bit-reversed;
 * @param r The number of bits to take into consideration when reversing.
 * @return The number k, bit-reversed according to integers with r bits.
 */
int bit_reverse(int k, int r);


/**
 * Smallest prime factor of a given number. If the argument is prime itself, then it is the return
 * value. This function is not exported, as it is only used internally.
 *
 * @param n Number to be inspected.
 * @return The smallest prime factor, or the number itself if it is already a prime.
 */
int factor(int n);


/**************************************************************************************************
 FFT implementations:
 **************************************************************************************************/
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
void direct_ft(Complex x[], Complex X[], int N);


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
void recursive_fft(Complex x[], Complex X[], int N);


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
void iterative_fft(Complex x[], Complex X[], int N);


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
void recursive_nfft(Complex x[], Complex X[], int N);


#endif
