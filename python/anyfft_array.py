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
# $ python3 anyfft_array.py


####################################################################################################
# Import needed modules:
from numpy import *                            # Deal with arrays;
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
    x = arange(size, dtype=complex)            # Generate a vector;
    t0 = perf_counter()                        # Start a timer;
    for j in range(0, repeat):                 # Repeated calls;
        f(x)
    return (perf_counter() - t0) / repeat      # Compute average;


####################################################################################################
# Direct FT:
def direct_ft(x):
    """
    Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    complexity. This implementation uses NumPy arrays for conciseness.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    X = zeros(x.shape, dtype=complex)          # Accumulate the results;
    W = exp(-2j*pi/N)                          # Twiddle factors;
    Wk = 1.
    for k in range(0, N):                      # Compute the kth coefficient;
        Wkn = 1.
        for n in range(0, N):                  #   Operate the summation;
            X[k] = X[k] + x[n]*Wkn             #     Compute every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        Wk = Wk * W
    return X


####################################################################################################
# Auxiliary function:
def __factor(n):
    """
    Smallest prime factor of a given number. If the argument is prime itself, then it is the return
    value.

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
    Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
    complexity than the direct FT, though the exact value is difficult to compute. This
    implementation uses NumPy arrays for conciseness.

    :Parameters:
      x
        The vector of which the FFT will be computed. It must be a composite number, or else the
        computation will be defered to the direct FT, and there will be no efficiency gain.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    N1 = __factor(N)                           # Find the smallest factor of the vector length;
    if N1 == N:                                # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else:
        N2 = N // N1                           # Decompose in two factors, N1 being prime;
        X = zeros((N, ), dtype=complex)        # Accumulate the results;
        W = exp(-2j*pi/N)                      # Twiddle factors;
        Wj = 1.
        for j in range(N1):                    # Compute every subsequence of size N2;
            Xj = recursive_fft(x[j::N1])       # Recursively compute the Fourier Transform;
            Wkj = 1.
            for k in range(N):
                X[k] = X[k] + Xj[k%N2] * Wkj   # Recombine results;
                Wkj = Wkj * Wj                 # Update twiddle factors;
            Wj = Wj * W
        return X


####################################################################################################
# Recursive FFT:
def vec_recursive_fft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
    complexity than the direct FT, though the exact value is difficult to compute. This
    implementation uses NumPy arrays for conciseness. In this implementation, loops are avoided by
    vectorizing the computation of the twiddle factors.

    :Parameters:
      x
        The vector of which the FFT will be computed. It must be a composite number, or else the
        computation will be defered to the direct FT, and there will be no efficiency gain.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    N1 = __factor(N)                           # Find the smallest factor of the vector length;
    if N1 == N:                                # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else:
        N2 = N // N1                           # Decompose in two factors, N1 being prime;
        X = zeros((N, ), dtype=complex)        # Accumulate the results;
        k = arange(N)
        for j in range(N1):                    # Compute every subsequence of size N2;
            Xj = vec_recursive_fft(x[j::N1])   # Recursively compute the Fourier Transform;
            Wkj = exp(-2j*pi*k*j/N)
            X = X + Xj[k%N2]*Wkj         # Recombine results;
        return X


####################################################################################################
# Main program:
if __name__ == "__main__":

    SIZES = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]

    # Start printing the table with time comparisons:
    print("+---------"*6 + "+")
    print("|    N    |   N^2   | Direct  | Recurs. | VecRec. | Intern. |")
    print("+---------"*6 + "+")

    # Try it with vectors with the given sizes:
    for n in SIZES:

        # Compute the average execution time:
        dtime  = time_it(direct_ft, n, REPEAT)
        rtime  = time_it(recursive_fft, n, REPEAT)
        vtime  = time_it(vec_recursive_fft, n, REPEAT)
        intime = time_it(fft.fft, n, REPEAT)

        # Print the results:
        tup = (n, n**2, dtime, rtime, intime)
        print(f'| {n:7} | {n**2:7} | {dtime:7.4f} | {rtime:7.4f} | {vtime:7.4f} | {intime:7.4f} |')

    print("+---------"*6 + "+")




