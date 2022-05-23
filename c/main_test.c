/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version compares Cooley-Tukey algorithm.
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o fft main_test.c fft.c fft_native.c my_complex.c time_it.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fft
 **************************************************************************************************/

// Include necessary libraries:
#include "test_it.h"                           // Time measurement;
#include "my_complex.h"                        // Mini-library for complex numbers;
#include "fft.h"                               // FFT implementations;
#include "fft_native.h"                        // FFT implementations using C native complex;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    Complex x[16], X[16];

    // Tests for the implementations of the FFT using `my_complex.h` library:
    printf("Direct FT - ");
    test_it(direct_ft, 8);
    printf("Recursive FFT - ");
    test_it(recursive_fft, 8);
    printf("Iterative FFT - ");
    test_it(iterative_fft, 8);
    printf("Direct FT - ");
    test_it(direct_ft, 16);
    printf("Recursive FFT - ");
    test_it(recursive_fft, 16);
    printf("Iterative FFT - ");
    test_it(iterative_fft, 16);

    // Tests for the implementations of the FFT using native complex:
    printf("Direct FT with Native Complex - ");
    native_complex_test_it(native_complex_direct_ft, 8);
    printf("Recursive FFT with Native Complex - ");
    native_complex_test_it(native_complex_recursive_fft, 8);
    printf("Iterative FFT with Native Complex - ");
    native_complex_test_it(native_complex_iterative_fft, 8);
    printf("Direct FT with Native Complex - ");
    native_complex_test_it(native_complex_direct_ft, 16);
    printf("Recursive FFT with Native Complex - ");
    native_complex_test_it(native_complex_recursive_fft, 16);
    printf("Iterative FFT with Native Complex - ");
    native_complex_test_it(native_complex_iterative_fft, 16);

    // Tests for the implementations of the FFT using `my_complex.h` library:
    printf("Direct FT - ");
    test_it(direct_ft, 12);
    printf("Recursive FFT - ");
    test_it(recursive_nfft, 12);

    // Tests for the implementations of the FFT using native complex:
    printf("Direct FT with Native Complex - ");
    native_complex_test_it(native_complex_direct_ft, 12);
    printf("Recursive FFT with Native Complex - ");
    native_complex_test_it(native_complex_recursive_nfft, 12);

    return 0;
}
