# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version
# This version compares Cooley-Tukey algorithm.
####################################################################################################
# Since Python is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ python main_test.py
#
# Notice that, since this implementation uses only standard Python 3 objects, you can try to run
# this script using Pypy 3, Cython 3 or any other Python 3 implementation.


####################################################################################################
# Import needed modules:
import fft_list
import fft_array
import fft_numpy
from test_it import test_it                    # Time events


####################################################################################################
# Main program:
if __name__ == "__main__":

    # Tests for the implementations of the FFT using Python native lists:
    print("Direct FT Using Lists - ")
    test_it(fft_list.direct_ft, 8)
    print("Direct FT Using Comprehension Lists - ")
    test_it(fft_list.lc_dft, 8)
    print("Recursive FFT Using Lists - ")
    test_it(fft_list.recursive_fft, 8)
    print("Iterative FFT Using Lists - ")
    test_it(fft_list.iterative_fft, 8)
    print("Direct FT Using Lists - ")
    test_it(fft_list.direct_ft, 16)
    print("Direct FT Using Comprehension Lists - ")
    test_it(fft_list.lc_dft, 16)
    print("Recursive FFT Using Lists - ")
    test_it(fft_list.recursive_fft, 16)
    print("Iterative FFT Using Lists - ")
    test_it(fft_list.iterative_fft, 16)

    # Tests for the implementations of the FFT using Python native arrays:
    print("Direct FT Using Native Arrays - ")
    test_it(fft_array.direct_ft, 8)
    print("Recursive FFT Using Native Arrays - ")
    test_it(fft_array.recursive_fft, 8)
    print("Iterative FFT Using Native Arrays - ")
    test_it(fft_array.iterative_fft, 8)
    print("Direct FT Using Native Arrays - ")
    test_it(fft_array.direct_ft, 16)
    print("Recursive FFT Using Native Arrays - ")
    test_it(fft_array.recursive_fft, 16)
    print("Iterative FFT Using Native Arrays - ")
    test_it(fft_array.iterative_fft, 16)

    # Tests for the implementations of the FFT using NumPy:
    print("Direct FT Using NumPy - ")
    test_it(fft_numpy.direct_ft, 8)
    print("Simplified FT Using NumPy - ")
    test_it(fft_numpy.simplified_ft, 8)
    print("Matrix FT Using NumPy - ")
    test_it(fft_numpy.matrix_ft, 8)
    print("Recursive FFT Using NumPy - ")
    test_it(fft_numpy.recursive_fft, 8)
    print("Iterative FFT Using NumPy - ")
    test_it(fft_numpy.iterative_fft, 8)
    print("Direct FT Using NumPy - ")
    test_it(fft_numpy.direct_ft, 16)
    print("Simplified FT Using NumPy - ")
    test_it(fft_numpy.simplified_ft, 16)
    print("Matrix FT Using NumPy - ")
    test_it(fft_numpy.matrix_ft, 16)
    print("Recursive FFT Using NumPy - ")
    test_it(fft_numpy.recursive_fft, 16)
    print("Iterative FFT Using NumPy - ")
    test_it(fft_numpy.iterative_fft, 16)

    # Tests for the implementations of the FFT using Python native lists:
    print("Direct FT Using Lists - ")
    test_it(fft_list.direct_ft, 12)
    print("Recursive FFT Using Lists - ")
    test_it(fft_list.recursive_nfft, 12)

    # Tests for the implementations of the FFT using Python native arrays:
    print("Direct FT Using Native Arrays - ")
    test_it(fft_array.direct_ft, 12)
    print("Recursive FFT Using Native Arrays - ")
    test_it(fft_array.recursive_nfft, 12)

    # Tests for the implementations of the FFT using NumPy:
    print("Direct FT Using NumPy - ")
    test_it(fft_list.direct_ft, 12)
    print("Recursive FFT Using NumPy - ")
    test_it(fft_list.recursive_nfft, 12)
