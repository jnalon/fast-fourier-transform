/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This version compares Cooley-Tukey algorithm for powers of 2 only.
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C++ compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ g++ -o fft main_fft.cpp fixed_point.cpp -lm
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

    // Start by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | N logN  | FDirect | FRecur. | FItera. ";
    cout << "| DDirect | DRecur. | DItera. | FPDir.  | FPRec.  | FPIte.  |" << endl;
    cout << "+---------+---------+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int r=5; r<11; r++) {

        // Compute the average execution time:
        int n = (int) exp2(r);
        float fdtime = time_it<float>(direct_ft<float>, n);
        float frtime = time_it<float>(recursive_fft<float>, n);
        float fitime = time_it<float>(iterative_fft<float>, n);
        float ddtime = time_it<double>(direct_ft<double>, n);
        float drtime = time_it<double>(recursive_fft<double>, n);
        float ditime = time_it<double>(iterative_fft<double>, n);
        float fpdtime = time_it<FixedPoint>(direct_ft<FixedPoint>, n);
        float fprtime = time_it<FixedPoint>(recursive_fft<FixedPoint>, n);
        float fpitime = time_it<FixedPoint>(iterative_fft<FixedPoint>, n);

        // Print the results:
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) <<   r*n << " ";
        cout << "| " << setw(7) << setprecision(4) <<  fdtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  frtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  fitime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  ddtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  drtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  ditime << " ";
        cout << "| " << setw(7) << setprecision(4) << fpdtime << " ";
        cout << "| " << setw(7) << setprecision(4) << fprtime << " ";
        cout << "| " << setw(7) << setprecision(4) << fpitime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+---------+---------+" << endl;
    return 0;
}
