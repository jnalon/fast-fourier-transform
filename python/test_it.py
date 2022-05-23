# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version
# This version implements a function to time repetitions of the FFT functions.
####################################################################################################

####################################################################################################
# Import needed modules:
from time import perf_counter                  # Time events;


####################################################################################################
# Definitions:
REPEATS = 50                                   # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function:
def test_it(f, size):
    """
    Pretty print of input and output of the Fourier transform for visual inspection.

    :Parameters:
      f
        Function to be called;
      size
        Number of elements in the vector on which the transform will be applied.
    """
    x = [ j+0j for j in range(size) ]
    X = f(x)
    print(f"N = {size} | Input | Output:")
    for i in range(size):
        print(f"  {i:2d}: ({x[i].real:8.4f}, {x[i].imag:8.4f}) | " +
              f"({X[i].real:8.4f}, {X[i].imag:8.4f})")
    print("-" * 30)


def time_it(f, size, repeats=REPEATS):
    """
    Measure execution time through repeated calls to a (Fast) Fourier Transform function.

    :Parameters:
      f
        Function to be called;
      size
        Size of the vector on which the transform will be applied;
      repeats
        Number of times the function will be called. Defaults to REPEAT.

    :Returns:
      The average execution time for that function with a vector of the given size.
    """
    x = [ j+0j for j in range(size) ]          # Generate a vector;
    t0 = perf_counter()                        # Start a timer;
    for j in range(0, repeats):                # Repeated calls;
        f(x)
    return (perf_counter() - t0) / repeats     # Compute average;
