/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This program performs a simple test of the implementations.
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C++ compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ g++ -o fft main_test.cpp fixed_point.cpp -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fft
 **************************************************************************************************/

// Include necessary libraries:
#include <iostream>                            // Input and Output;
#include "my_complex.h"
#include "fixed_point.h"
#include "test_it.h"
#include "fft.h"


using namespace std;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    // Tests for the implementation using complex floats:
    cout << "Direct FT with Float Complex - ";
    test_it<float>(direct_ft<float>, 8);
    cout << "Recursive FFT with Float Complex - ";
    test_it<float>(recursive_fft<float>, 8);
    cout << "Iterative FFT with Float Complex - ";
    test_it<float>(iterative_fft<float>, 8);
    cout << "Direct FT with Float Complex - ";
    test_it<float>(direct_ft<float>, 16);
    cout << "Recursive FFT with Float Complex - ";
    test_it<float>(recursive_fft<float>, 16);
    cout << "Iterative FFT with Float Complex - ";
    test_it<float>(iterative_fft<float>, 16);

    // Tests for the implementation using complex doubles:
    cout << "Direct FT with Double Complex - ";
    test_it<double>(direct_ft<double>, 8);
    cout << "Recursive FFT with Double Complex - ";
    test_it<double>(recursive_fft<double>, 8);
    cout << "Iterative FFT with Double Complex - ";
    test_it<double>(iterative_fft<double>, 8);
    cout << "Direct FT with Double Complex - ";
    test_it<double>(direct_ft<double>, 16);
    cout << "Recursive FFT with Double Complex - ";
    test_it<double>(recursive_fft<double>, 16);
    cout << "Iterative FFT with Double Complex - ";
    test_it<double>(iterative_fft<double>, 16);

    // Tests for the implementation using FixedPoint complex:
    cout << "Direct FT with FixedPoint Complex - ";
    test_it<FixedPoint>(direct_ft<FixedPoint>, 8);
    cout << "Recursive FFT with FixedPoint Complex - ";
    test_it<FixedPoint>(recursive_fft<FixedPoint>, 8);
    cout << "Iterative FFT with FixedPoint Complex - ";
    test_it<FixedPoint>(iterative_fft<FixedPoint>, 8);
    cout << "Direct FT with FixedPoint Complex - ";
    test_it<FixedPoint>(direct_ft<FixedPoint>, 16);
    cout << "Recursive FFT with FixedPoint Complex - ";
    test_it<FixedPoint>(recursive_fft<FixedPoint>, 16);
    cout << "Iterative FFT with FixedPoint Complex - ";
    test_it<FixedPoint>(iterative_fft<FixedPoint>, 16);

    // Tests for the implementation using complex floats:
    cout << "Direct FT with Float Complex - ";
    test_it<float>(direct_ft<float>, 12);
    cout << "Recursive FFT with Float Complex - ";
    test_it<float>(recursive_nfft<float>, 12);

    // Tests for the implementation using complex doubles:
    cout << "Direct FT with Double Complex - ";
    test_it<double>(direct_ft<double>, 12);
    cout << "Recursive FFT with Double Complex - ";
    test_it<double>(recursive_nfft<double>, 12);

    // Tests for the implementation using complex doubles:
    cout << "Direct FT with FixedPoint Complex - ";
    test_it<FixedPoint>(direct_ft<FixedPoint>, 12);
    cout << "Recursive FFT with FixedPoint Complex - ";
    test_it<FixedPoint>(recursive_nfft<FixedPoint>, 12);

}
