# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
####################################################################################################
# Since Python is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ python main_fft_pure.py
#
# Notice that, since this implementation uses only standard Python 3 objects, you can try to run
# this script using Pypy 3, Cython 3 or any other Python 3 implementation.


####################################################################################################
# Import needed modules:
import fft_list
import fft_array
from test_it import time_it                    # Time events;


####################################################################################################
# Main program:
if __name__ == "__main__":

    # Start by printing the table with time comparisons:
    print("+---------"*10 + "+")
    print("|    N    |   N^2   | N logN  "
          "| Direct  | CList   | Array   "
          "| Recurs. | ARecur. "
          "| Itera.  | AItera  |")
    print("+---------"*10 + "+")

    # Try it with vectors with size ranging from 32 to 1024 samples:
    for r in range(5, 11):

        # Compute the average execution time:
        n = 2**r
        dtime  = time_it(fft_list.direct_ft, n)
        ctime  = time_it(fft_list.lc_dft, n)
        atime  = time_it(fft_array.direct_ft, n)
        rtime  = time_it(fft_list.recursive_fft, n)
        artime = time_it(fft_array.recursive_fft, n)
        itime  = time_it(fft_list.iterative_fft, n)
        aitime = time_it(fft_array.iterative_fft, n)

        # Print the results:
        print(f'| {n:7} | {n**2:7} | {r*n:7} '
              f'| {dtime:7.4f} | {ctime:7.4f} | {atime:7.4f} '
              f'| {rtime:7.4f} | {artime:7.4f} '
              f'| {itime:7.4f} | {aitime:7.4f} |')

    print("+---------"*10 + "+")
