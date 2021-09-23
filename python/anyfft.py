# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version
# This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
#
# José Alexandre Nalon
####################################################################################################
# Since Python is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ python anyfft_pure.py
#
# Notice that, since this implementation uses only standard Python 3 objects, you can try to run
# this script using Pypy 3, Cython 3 or any other Python 3 implementation.


####################################################################################################
# Import needed modules:
import numpy as np                             # Deal with arrays;
import fft_list
import fft_array
import fft_numpy
from time import perf_counter                  # Time events;


####################################################################################################
# Definitions:
REPEAT = 50                                    # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function:
def time_it(f, size, repeat=REPEAT):
    """
    Measure execution time through repeated calls to a (Fast) Fourier Transform function.

    :Parameters:
      f
        Function to be called;
      size
        Size of the vector on which the transform will be applied;
      repeat
        Number of times the function will be called. Defaults to REPEAT.

    :Returns:
      The average execution time for that function with a vector of the given size.
    """
    x = [ j+0j for j in range(size) ]          # Generate a vector;
    t0 = perf_counter()                        # Start a timer;
    for j in range(0, repeat):                 # Repeated calls;
        f(x)
    return (perf_counter() - t0) / repeat      # Compute average;


####################################################################################################
# Main program:
if __name__ == "__main__":

    SIZES = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]

    # Start printing the table with time comparisons:
    print("+---------"*7 + "+")
    print("|    N    |   N^2   | Direct  | ADirect | Recurs. | ARecur. | RNumpy  |")
    print("+---------"*7 + "+")

    # Try it with vectors with the given sizes:
    for n in SIZES:

        # Compute the average execution time:
        dtime  = time_it(fft_list.direct_ft, n, REPEAT)
        atime  = time_it(fft_array.direct_ft, n, REPEAT)
        rtime  = time_it(fft_list.recursive_nfft, n, REPEAT)
        artime = time_it(fft_array.recursive_nfft, n, REPEAT)
        ntime  = time_it(fft_numpy.recursive_nfft, n, REPEAT)

        # Print the results:
        print(f'| {n:7} | {n**2:7} | {dtime:7.4f} |'
              f' {atime:7.4f} | {rtime:7.4f} | {artime:7.4f} | {ntime:7.4f} |')

    print("+---------"*7 + "+")


