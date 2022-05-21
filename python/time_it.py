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
