/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version compares Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o anyfft main_anyfft.c fft.c fft_float.c fft_double.c my_complex.c time_it.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./anyfft
 **************************************************************************************************/

// Include necessary libraries:
#include "test_it.h"                           // Time measurement;
#include "fft.h"                               // FFT implementations;
#include "fft_float.h"                         // FFT implementations using C float complex;
#include "fft_double.h"                        // FFT implementations using C double complex;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, fdtime, frtime, ddtime, drtime;
    int SIZES[] = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
    int n;

    // Start by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | Direct  | Recurs. | FDirect | FRecur. | DDirect | DRecur. |\n");
    printf("+---------+---------+---------+---------+---------+---------+---------+---------+\n");

    // Try it with vectors with the given sizes:
    for(int i=0; i<8; i++) {

        // Compute the average execution time:
        n = SIZES[i];
        dtime = time_it(direct_ft, n, REPEATS);
        rtime = time_it(recursive_nfft, n, REPEATS);
        fdtime = float_time_it(float_direct_ft, n, REPEATS);
        frtime = float_time_it(float_recursive_nfft, n, REPEATS);
        ddtime = double_time_it(double_direct_ft, n, REPEATS);
        drtime = double_time_it(double_recursive_nfft, n, REPEATS);

        // Print the results:
        printf("| %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, dtime, rtime, fdtime, frtime, ddtime, drtime);
    }
    printf("+---------+---------+---------+---------+---------+---------+---------+---------+\n");

    return 0;
}

