# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version using native lists.
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# This is a function only library, and executing this file won't give any results. Notice however,
# that, since this implementation uses only standard Python 3 objects, you can try to include this
# in scripts and run them using Pypy 3, Cython 3 or any other Python 3 implementation.


####################################################################################################
# Import needed modules:
import math                                    # Mathematical operations;
import cmath as cm                             # Complex math;


####################################################################################################
# Direct FT:
def direct_ft(x):
    """
    Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    complexity. This implementation uses native Python lists.

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
    W = cm.exp(-2j*cm.pi/N)                    # Twiddle factors;
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
lc_dft = lambda x: [ sum([ xn*cm.exp(-2j*cm.pi*n*k/len(x)) for n, xn in enumerate(x) ])
                     for k, _ in enumerate(x) ]


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
        W = [ cm.exp(-2j*cm.pi*k/N) for k in range(N//2) ]     # Twiddle factors;
        WXo = [ Wk*Xok for Wk, Xok in zip(W, Xo) ]
        X = ([ Xek + WXok for Xek, WXok in zip(Xe, WXo) ] +    # Recombine results;
             [ Xek - WXok for Xek, WXok in zip(Xe, WXo) ])
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
    r = int(math.log2(N))                      # Number of bits;
    X = [ complex(xi) for xi in x ]            # Accumulate the results;
    for k in range(0, N):
        l = bit_reverse(k, r)                  # Reorder the vector according to the
        X[l] = x[k]                            #   bit-reversed order;

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k in range(0, r):
        for l in range(0, N, 2*step):
            W = cm.exp(-1j * cm.pi / step)     # Twiddle factors;
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
    rn = int(math.ceil(cm.sqrt(n).real))       # Search up to the square root of the number;
    for i in range(2, rn+1):
        if n%i == 0:                           # When remainder is zero, factor is found;
            return i
    return n


####################################################################################################
# Recursive FFT:
def recursive_nfft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm for vectors of length
    different from a power of 2. This has smaller complexity than the direct FT, though the exact
    value is difficult to compute. This implementation uses native Python lists.

    :Parameters:
      x
        The vector of which the FFT will be computed. Its length must be a composite number, or else
        the computation will be defered to the direct FT, and there will be no efficiency gain.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    N1 = __factor(N)                           # Find the smallest factor of the vector length;
    if N1 == N:                                # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else:
        N2 = N // N1                           # Decompose in two factors, N1 being prime;
        X = [ 0+0j ] * N                       # Accumulate the results;
        W = cm.exp(-2j*cm.pi/N)                # Twiddle factors;
        Wj = 1.
        for j in range(N1):                    # Compute every subsequence of size N2;
            Xj = recursive_nfft(x[j::N1])
            Wkj = 1.
            for k in range(N):
                X[k] = X[k] + Xj[k%N2] * Wkj   # Recombine results;
                Wkj = Wkj * Wj                 # Update twiddle factors;
            Wj = Wj * W
        return X


####################################################################################################
# Main program:
if __name__ == "__main__":
    pass
