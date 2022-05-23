/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the header file for the time measuring functions.
 **************************************************************************************************/

#ifndef __TIME_IT__
#define __TIME_IT__


// Include necessary libraries:
#include <stdlib.h>                            // Standard Library;
#include <complex.h>                           // Complex numbers;
#include <time.h>                              // Time measurement;
#include "my_complex.h"                        // Mini-library for complex numbers;


// Definitions:
#define REPEATS 500


/**
 * Pretty print of input and output of the Fourier transform for visual inspection.
 *
 * @param f Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 * @param size Number of elements in the vector on which the transform will be applied;
 */
void test_it(void (*f)(Complex *, Complex *, int), int size);


/**
 * Pretty print of input and output of the Fourier transform for visual inspection, using native
 * complex number library.
 *
 * @param f Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 * @param size Number of elements in the vector on which the transform will be applied;
 */
void native_complex_test_it(void (*f)(float complex *, float complex *, int), int size);


/**
 * Measure execution time through repeated calls to a (Fast) Fourier Transform function.
 *
 * @param f Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 * @param size Number of elements in the vector on which the transform will be applied;
 * @param repeat Number of times the function will be called.
 * @return The average execution time for that function with a vector of the given size.
 */
float time_it(void (*f)(Complex *, Complex *, int), int size, int repeat);


/**
 * Measure execution time through repeated calls to a (Fast) Fourier Transform function for
 * functions using native C complex numbers.
 *
 * @param f Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 * @param size Number of elements in the vector on which the transform will be applied;
 * @param repeat Number of times the function will be called.
 * @return The average execution time for that function with a vector of the given size.
 */
float native_complex_time_it(void (*f)(float complex *, float complex *, int), int size, int repeat);


#endif
