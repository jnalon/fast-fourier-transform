# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# Since Python is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ python fft.py


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
        Power of two of the size of the vector on which the transform will be applied;
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
    complexity. This implementation uses native Python lists -- apparently, it does have an impact
    on code, resulting in faster execution.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    X = [ 0+0j ] * N                           # Accumulate the results;
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
# This is a small and very readable DFT implementation using comprehension lists. This shows how
# much of Python can be directly translated from Math directly to code. An explanation of every term
# follows:
#
# lambda x:                                    # Given in functional form;
#   [ sum(                                     # Summation
#       [
#         xn                                   #   of every sample
#         * exp(-2j*pi*n*k/len(x))             #   times the twiddle factor
#             for n, xn in enumerate(x)        #   over the interval of samples;
#       ])
#       for k, _ in enumerate(x) ]             # Repeat for every coefficient;
#
# Here is the analysis equation (in LaTeX form) for comparison:
#
#  X[k] = \sum_{0}^{N-1} x[n] e^{-j 2 \pi k n / N} \; k = 0 \ldots N-1
#
# This implementation is not, of course, very efficient, since it doesn't take advantage of the
# regularity of the twiddle factors, and computes complex exponentials for every term.
lc_dft = lambda x: [ sum([ xn*exp(-2j*pi*n*k/len(x)) for n, xn in enumerate(x) ]) for k, _ in enumerate(x) ]


####################################################################################################
# Direct FT using NumPy arrays:
def array_direct_ft(x):
    """
    Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    complexity. This implementation uses NumPy arrays for memory conservation.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

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
# Recursive Decimation-in-time FFT:
def recursive_fft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
    complexity. This implementation uses native Python lists.

    :Parameters:
      x
        The vector of which the FFT will be computed. This should always be called with a vector of
        a power of two length, or it will fail. No checks on this are made.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    if len(x) == 1:                            # A length-1 vector is its own FT;
        return x
    else:
        N = len(x)                             # Length of the vector;
        Xe = recursive_fft(x[0::2])            # Transform of even samples;
        Xo = recursive_fft(x[1::2])            # Transform of odd samples;
        W = [ exp(-2j*pi*k/N) for k in range(0, int(N/2)) ]    # Twiddle factors;
        WXo = [ Wk*Xok for Wk, Xok in zip(W, Xo) ]
        X = [ Xek + WXok for Xek, WXok in zip(Xe, WXo) ] + \
            [ Xek - WXok for Xek, WXok in zip(Xe, WXo) ]       # Recombine results;
        return X


####################################################################################################
# Recursive Decimation-in-time FFT using NumPy arrays:
def array_recursive_fft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
    complexity. This implementation uses NumPy arrays.

    :Parameters:
      x
        The vector of which the FFT will be computed. This should always be called with a vector of
        a power of two length, or it will fail. No checks on this are made.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    if len(x) == 1:                                    # A length-1 vector is its own FT;
        return x
    else:
        N = len(x)                                     # Length of the vector;
        Xe = array_recursive_fft(x[0::2])              # Transform of even samples;
        Xo = array_recursive_fft(x[1::2])              # Transform of odd samples;
        W = exp(-2j*pi/N) ** range(0, int(N/2))        # Twiddle factors;
        WXo = W * Xo                                   # Repeated computation;
        X = hstack((Xe + WXo, Xe - WXo))               # Recombine results;
        return X


####################################################################################################
# Auxililary function to reorder a vector in bit-reversed order:
def bit_reverse(k, r):
    """
    Bit-reversed version of an integer number.

    :Parameters:
      k
        The number to be bit-reversed;
      r
        The number of bits to take into consideration when reversing.

    :Returns:
      The number k, bit-reversed according to integers with r bits.
    """
    l = 0                                      # Accumulate the results;
    for i in range(0, r):                      # Loop on every bit;
        l = (l << 1) + (k & 1)                 # Test less signficant bit and add;
        k = (k >> 1)                           # Test next bit;
    return l


####################################################################################################
# Iterative Decimation-in-time FFT, using lists:
def iterative_fft(x):
    """
    Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
    O(N log_2(N)) complexity, and since there are less function calls, it will probably be
    marginally faster than the recursive versions. It uses native Python lists.

    :Parameters:
      x
        The vector of which the FFT will be computed. This should always be called with a vector of
        a power of two length, or it will fail. No checks on this are made.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of vector;
    r = int(log(N)/log(2))                     # Number of bits;
    X = [ complex(xi) for xi in x ]            # Accumulate the results;
    for k in range(0, N):
        l = bit_reverse(k, r)                  # Reorder the vector according to the
        X[l] = x[k]                            #   bit-reversed order;

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k in range(0, r):
        for l in range(0, N, 2*step):
            W = exp(-1j * pi / step)           # Twiddle factors;
            Wkn = 1.
            for n in range(0, step):
                p = l + n
                q = p + step
                X[q] = X[p] - Wkn*X[q]         # Recombine results;
                X[p] = 2*X[p] - X[q]
                Wkn = Wkn * W                  # Update twiddle factors;
        step = 2*step

    return X


####################################################################################################
# Iterative Decimation-in-time FFT, using NumPy arrays:
def array_iterative_fft(x):
    """
    Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
    O(N log_2(N)) complexity, and since there are less function calls, it will probably be
    marginally faster than the recursive versions. It uses NumPy arrays.

    :Parameters:
      x
        The vector of which the FFT will be computed. This should always be called with a vector of
        a power of two length, or it will fail. No checks on this are made.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of vector;
    r = int(log(N)/log(2))                     # Number of bits;
    X = array(x, dtype=complex)                # Accumulate the results;
    for k in range(0, N):
        l = bit_reverse(k, r)                  # Reorder the vector according to the
        X[l] = x[k]                            #   bit-reversed order;

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k in range(0, r):
        for l in range(0, N, 2*step):
            W = exp(-1j * pi / step)           # Twiddle factors;
            Wkn = 1.
            for n in range(0, step):
                p = l + n
                q = p + step
                X[q] = X[p] - Wkn*X[q]         # Recombine results;
                X[p] = 2*X[p] - X[q]
                Wkn = Wkn * W                  # Update twiddle factors;
        step = 2*step

    return X


####################################################################################################
# Main program:
if __name__ == "__main__":

    # Start by printing the table with time comparisons:
    print("+---------"*11 + "+")
    print("|    N    |   N^2   | N logN  | Direct  | CList   | Array   | Recurs. | Rarray  | Itera.  | AItera  | Interna |")
    print("+---------"*11 + "+")

    # Try it with vectors with size ranging from 32 to 1024 samples:
    for r in range(5, 11):

        # Compute the average execution time:
        n = 2**r
        dtime  = time_it(direct_ft, n, REPEAT)
        ctime  = time_it(lc_dft, n, REPEAT)
        atime  = time_it(array_direct_ft, n, REPEAT)
        rtime  = time_it(recursive_fft, n, REPEAT)
        artime = time_it(array_recursive_fft, n, REPEAT)
        itime  = time_it(iterative_fft, n, REPEAT)
        aitime = time_it(array_iterative_fft, n, REPEAT)
        intime = time_it(fft.fft, n, REPEAT)

        # Print the results:
        tup = (n, n**2, r*n, dtime, ctime, atime, rtime, artime, itime, aitime, intime)
        print(f'| {n:7} | {n**2:7} | {r*n:7} | {dtime:7.4f} | {ctime:7.4f} | {atime:7.4f} ' \
              f'| {rtime:7.4f} | {artime:7.4f} | {itime:7.4f} | {aitime:7.4f} | {intime:7.4f} |')

    print("+---------"*11 + "+")
