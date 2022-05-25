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

    vector<int> sizes = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };

    // Start by printing the table with time comparisons:
    cout << "+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+" << endl;
    cout << "|    N    |   N^2   | FDirect | FRecur. ";
    cout << "| DDirect | DRecur. | FPDir.  | FPRec.  |" << endl;
    cout << "+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+" << endl;

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for(int n : sizes) {

        // Compute the average execution time:
        float fdtime = time_it<float>(direct_ft<float>, n);
        float frtime = time_it<float>(recursive_nfft<float>, n);
        float ddtime = time_it<double>(direct_ft<double>, n);
        float drtime = time_it<double>(recursive_nfft<double>, n);
        float fpdtime = time_it<FixedPoint>(direct_ft<FixedPoint>, n);
        float fprtime = time_it<FixedPoint>(recursive_nfft<FixedPoint>, n);

        // Print the results:
        cout << "| " << setw(7) <<     n << " ";
        cout << "| " << setw(7) <<   n*n << " ";
        cout << "| " << setw(7) << setprecision(4) <<  fdtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  frtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  ddtime << " ";
        cout << "| " << setw(7) << setprecision(4) <<  drtime << " ";
        cout << "| " << setw(7) << setprecision(4) << fpdtime << " ";
        cout << "| " << setw(7) << setprecision(4) << fprtime << " |" << endl;
    }

    cout << "+---------+---------+---------+---------";
    cout << "+---------+---------+---------+---------+" << endl;
    return 0;
}
