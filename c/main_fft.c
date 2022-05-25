/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This version compares Cooley-Tukey algorithm for powers of 2 only.
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ gcc -o fft main_fft.c fft.c fft_float.c fft_double.c my_complex.c time_it.c -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fft
 **************************************************************************************************/

// Include necessary libraries:
#include "test_it.h"                           // Time measurement;
#include "my_complex.h"                        // Mini-library for complex numbers;
#include "fft.h"                               // FFT implementations;
#include "fft_float.h"                         // FFT implementations using native C float complex;
#include "fft_double.h"                        // FFT implementations using native C double complex;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    float dtime, rtime, itime, fdtime, frtime, fitime, ddtime, drtime, ditime;
    int n;

    // Start by printing the table with time comparisons:
    printf("+---------+---------+---------+---------+---------+---------+"
            "---------+---------+---------+---------+---------+---------+\n");
    printf("|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  |"
            " FDirect | FRecur. | FItera. | DDirect | DRecur. | DItera. |\n");
    printf("+---------+---------+---------+---------+---------+---------+"
            "---------+---------+---------+---------+---------+---------+\n");

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int r=5; r<11; r++) {

        // Compute the average execution time:
        n = exp2(r);
        dtime = time_it(direct_ft, n, REPEATS);
        rtime = time_it(recursive_fft, n, REPEATS);
        itime = time_it(iterative_fft, n, REPEATS);
        fdtime = float_time_it(float_direct_ft, n, REPEATS);
        frtime = float_time_it(float_recursive_fft, n, REPEATS);
        fitime = float_time_it(float_iterative_fft, n, REPEATS);
        ddtime = double_time_it(double_direct_ft, n, REPEATS);
        drtime = double_time_it(double_recursive_fft, n, REPEATS);
        ditime = double_time_it(double_iterative_fft, n, REPEATS);

        // Print the results:
        printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, r*n, dtime, rtime, itime, fdtime, frtime, fitime, ddtime, drtime, ditime);
    }
    printf("+---------+---------+---------+---------+---------+---------+"
            "---------+---------+---------+---------+---------+---------+\n");

    return 0;
}
