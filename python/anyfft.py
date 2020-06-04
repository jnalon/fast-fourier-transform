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
# $ python anyfft.py


####################################################################################################
# Import needed modules:
from numpy import *                            # Deals with arrays;
from time import perf_counter                  # Time events;


####################################################################################################
# Definitions:
REPEAT = 50                                    # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function:
def time_it(f, size, repeat=REPEAT):
    """
    This function calls a Fast Fourier Transform function repeatedly a certain number of times,
    measure execution time and average it.

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
    x = arange(size, dtype=complex)            # Generate a vector;
    t0 = perf_counter()                        # Starts a timer;
    for j in range(0, repeat):                 # Repeated calls;
        f(x)
    return (perf_counter() - t0) / repeat      # Computes average;


####################################################################################################
# Direct FT:
def direct_ft(x):
    """
    Computes the Discrete Fourier Transform directly from the definition, an algorithm that has
    O(N^2) complexity. This implementation uses native Python lists -- apparently, it does have an
    impact on code, resulting in faster execution.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    X = [ 0+0j ] * N                           # Accumulates the results;
    W = exp(-2j*pi/N)                          # Twiddle factors;
    Wk = 1.
    for k in range(0, N):                      # Compute the kth coefficient;
        Wkn = 1.
        for n in range(0, N):                  #   Operates the summation;
            X[k] = X[k] + x[n]*Wkn             #     Computes every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        Wk = Wk * W
    return X


####################################################################################################
# Auxiliary function:
def __factor(n):
    """
    This function finds the smallest prime factor of a given number. If the argument is prime
    itself, then it is the return value.

    :Parameters:
      n
        Number to be inspected.

    :Returns:
      The smallest prime factor, or the number itself if it is already a prime.
    """
    rn = int(ceil(sqrt(n)))                    # Search up to the square root of the number;
    for i in range(2, rn+1):
        if n%i == 0:                           # When remainder is zero, factor is found;
            return i
    return n


####################################################################################################
# Recursive FFT:
def recursive_fft(x):
    """
    Computes the Fast Fourier Transform using a recursive decimation in time algorithm. This has
    smallest complexity than the direct FT, though the exact value is difficult to compute. This
    implementation uses NumPy arrays for conciseness.

    :Parameters:
      x
        The vector of which the FFT will be computed. It must be a composite number, or else the
        computation will be defered to the direct FT, and there will be no efficiency gain.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    N1 = __factor(N)                           # Finds the smallest factor of the vector length;
    if N1 == N:                                # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else:
        N2 = N // N1                           # Decompose in two factors, N1 being prime;
        X = zeros((N, ), dtype=complex)        # Accumulates the results;
        W = exp(-2j*pi/N)                      # Twiddle factors;
        Wj = 1.
        for j in range(N1):                    # Computes every subsequence of size N2;
            Xj = recursive_fft(x[j::N1])
            Wkj = 1.
            for k in range(N):                 # Recombine results;
                X[k] = X[k] + Xj[k%N2] * Wkj
                Wkj = Wkj * Wj                 # Update twiddle factors;
            Wj = Wj * W
        return X


####################################################################################################
# Main program:
if __name__ == "__main__":

    # Starts printing the table with time comparisons:
    print("+---------"*5 + "+")
    print("|    N    |   N^2   | Direct  | Recurs. | Interna |")
    print("+---------"*5 + "+")

    # Try it with vectors with the given sizes:
    sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]
    for n in sizes:

        # Computes the average execution time:
        dtime  = time_it(direct_ft, n, REPEAT)
        rtime  = time_it(recursive_fft, n, REPEAT)
        intime = time_it(fft.fft, n, REPEAT)

        # Print the results:
        tup = (n, n**2, dtime, rtime, intime)
        print(f'| {n:7} | {n**2:7} | {dtime:7.4f} | {rtime:7.4f} | {intime:7.4f} |')

    print("+---------"*5 + "+")
