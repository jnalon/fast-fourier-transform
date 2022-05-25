/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This is the implementation file for the testing and time measuring functions.
 **************************************************************************************************/

#ifndef __TIME_IT__
#define __TIME_IT__


// Include necessary libraries:
#include <iostream>                            // Input and Output;
#include <iomanip>                             // I/O Manipulation;
#include <vector>                              // Deals with arrays;
#include <chrono>                              // Time measurement;


using namespace std;


// Definitions:
const int repeats = 500;                       // Number of executions to compute average time;


/**
 * Pretty print of input and output of the Fourier transform for visual inspection.
 *
 * @param f Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 * @param size Number of elements in the vector on which the transform will be applied;
 */
template<typename T>
void test_it(vector<Complex<T>> (*f)(vector<Complex<T>>), int size)
{
    vector<Complex<T>> x(size), X;

    for(int i=0; i<size; i++)
        x[i] = Complex<T>(i, 0);
    X = f(x);
    std::cout << "N = " << size << " | Input | Output:" << std::endl;
    for(int i=0; i<size; i++)
        std::cout << "  " << i << ": " << x[i] << " | " << X[i] << std::endl;
}


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
template<typename T>
float time_it(vector<Complex<T>> (*f)(vector<Complex<T>>), int size, int repeat=repeats)
{
    vector<Complex<T>> x(size), X;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = Complex<T>(j, 0);
    auto t0 = chrono::steady_clock::now();     // Start counting time;
    for(int j=0; j<repeat; j++)
        X = (*f)(x);
    auto t1 = chrono::steady_clock::now();     // End of time measuring;
    return chrono::duration_cast<chrono::seconds>(t1 - t0).count() / (float) repeat;
}


#endif
