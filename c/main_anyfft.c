/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 * José Alexandre Nalon
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o anyfft main_anyfft.c fft.c fft_native.c my_complex.c time_it.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./anyfft
 **************************************************************************************************/

// Include necessary libraries:
#include "time_it.h"                           // Time measurement;
#include "fft.h"                               // FFT implementations;
#include "fft_native.h"                        // FFT implementations using C native complex;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, ncdtime, ncrtime;
    int SIZES[] = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };
    int n;

    // Start by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | Direct  | Recurs. | NDirect | NRecur. |\n");
    printf("+---------+---------+---------+---------+---------+---------+\n");

    // Try it with vectors with the given sizes:
    for(int i=0; i<8; i++) {

        // Compute the average execution time:
        n = SIZES[i];
        dtime = time_it(direct_ft, n, REPEATS);
        rtime = time_it(recursive_nfft, n, REPEATS);
        ncdtime = native_complex_time_it(native_complex_direct_ft, n, REPEATS);
        ncrtime = native_complex_time_it(native_complex_recursive_nfft, n, REPEATS);

        // Print the results:
        printf("| %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, dtime, rtime, ncdtime, ncrtime);
    }
    printf("+---------+---------+---------+---------+---------+---------+\n");

    return 0;
}

