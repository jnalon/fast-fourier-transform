# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version using native arrays.
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# This is a function only library, and executing this file won't give any results. Notice however,
# that, since this implementation uses only standard Python 3 objects, you can try to include this
# in scripts and run them using Pypy 3, Cython 3 or any other Python 3 implementation.


####################################################################################################
# Import needed modules:
from array import array                        # Efficient arrays;
import math                                    # Math operations;


####################################################################################################
# Direct FT:
def direct_ft(x):
    """
    Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    complexity. This implementation uses arrays from the standard Python library.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of the vector;
    x_real = array('f', [ xr.real for xr in x ])
    x_imag = array('f', [ xi.imag for xi in x ])
    X_real = array('f', [ 0. ] * N)            # Accumulate the results;
    X_imag = array('f', [ 0. ] * N)
    W_real = math.cos(-2*math.pi/N)            # Twiddle factors;
    W_imag = math.sin(-2*math.pi/N)
    Wk_real = 1.
    Wk_imag = 0.
    for k in range(0, N):                      # Compute the kth coefficient;
        Wkn_real = 1.
        Wkn_imag = 0.
        for n in range(0, N):                  #   Operate the summation;
            X_real[k] += x_real[n]*Wkn_real - x_imag[n]*Wkn_imag       # Compute every term;
            X_imag[k] += x_real[n]*Wkn_imag + x_imag[n]*Wkn_real
            aux_real = Wkn_real * Wk_real - Wkn_imag * Wk_imag         # Update twiddle factors;
            aux_imag = Wkn_real * Wk_imag + Wkn_imag * Wk_real
            Wkn_real = aux_real
            Wkn_imag = aux_imag
        aux_real = Wk_real * W_real - Wk_imag * W_imag
        aux_imag = Wk_real * W_imag + Wk_imag * W_real
        Wk_real = aux_real
        Wk_imag = aux_imag
    return [ Xr + 1j*Xi for Xr, Xi in zip(X_real, X_imag) ]


####################################################################################################
# Recursive Decimation-in-time FFT:
def recursive_fft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
    complexity. This implementation uses arrays from the standard Python library.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    if len(x) == 1:                            # A length-1 vector is its own FT;
        return x
    else:
        N = len(x)                             # Length of the vector;

        # Separated transform of even and odd samples:
        Xe = recursive_fft(x[0::2])
        Xo = recursive_fft(x[1::2])
        Xe_real = array('f', [ xr.real for xr in Xe ])
        Xe_imag = array('f', [ xi.imag for xi in Xe ])
        Xo_real = array('f', [ xr.real for xr in Xo ])
        Xo_imag = array('f', [ xi.imag for xi in Xo ])

        # Twiddle factors:
        kn = list(range(0, N//2))
        W_real = array('f', [ math.cos(-2*math.pi*k/N) for k in kn ])
        W_imag = array('f', [ math.sin(-2*math.pi*k/N) for k in kn ])

        # Multiply transform of odd samples by twiddle factors:
        tXo = list(zip(W_real, W_imag, Xo_real, Xo_imag))
        WXo_real = array('f', [ Wkr*Xokr - Wki*Xoki for Wkr, Wki, Xokr, Xoki in tXo ])
        WXo_imag = array('f', [ Wkr*Xoki + Wki*Xokr for Wkr, Wki, Xokr, Xoki in tXo ])

        # Aggregate results:
        a_real = list(zip(Xe_real, WXo_real))
        a_imag = list(zip(Xe_imag, WXo_imag))
        X_real = ( array('f', [ Xekr + WXokr for Xekr, WXokr in a_real ]) +
                   array('f', [ Xekr - WXokr for Xekr, WXokr in a_real ]) )
        X_imag = ( array('f', [ Xeki + WXoki for Xeki, WXoki in a_imag ]) +
                   array('f', [ Xeki - WXoki for Xeki, WXoki in a_imag ]) )
        return [ Xr + 1j*Xi for Xr, Xi in zip(X_real, X_imag) ]


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
    marginally faster than the recursive versions. This uses arrays from the standard Python
    library.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    N = len(x)                                 # Length of vector;
    r = int(math.log2(N))                      # Number of bits;
    X_real = array('f', [ xr.real for xr in x ])
    X_imag = array('f', [ xi.imag for xi in x ])
    for k in range(0, N):
        l = bit_reverse(k, r)                  # Reorder the vector according to the
        X_real[l] = x[k].real                  #   bit-reversed order;
        X_imag[l] = x[k].imag

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k in range(0, r):
        for l in range(0, N, 2*step):
            W_real = math.cos(-math.pi / step)         # Twiddle factors;
            W_imag = math.sin(-math.pi / step)
            Wkn_real = 1.
            Wkn_imag = 0.
            for n in range(0, step):
                p = l + n
                q = p + step

                # Recombine results:
                aux_real = Wkn_real*X_real[q] - Wkn_imag*X_imag[q]
                aux_imag = Wkn_real*X_imag[q] + Wkn_imag*X_real[q]
                X_real[q] = X_real[p] - aux_real
                X_imag[q] = X_imag[p] - aux_imag

                aux_real = 2*X_real[p] - X_real[q]
                aux_imag = 2*X_imag[p] - X_imag[q]
                X_real[p] = aux_real
                X_imag[p] = aux_imag

                # Update twiddle factors:
                aux_real = Wkn_real * W_real - Wkn_imag * W_imag
                aux_imag = Wkn_real * W_imag + Wkn_imag * W_real
                Wkn_real = aux_real
                Wkn_imag = aux_imag

        step = 2*step

    return [ Xr + 1j*Xi for Xr, Xi in zip(X_real, X_imag) ]


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
    rn = int(math.ceil(math.sqrt(n)))          # Search up to the square root of the number;
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
    value is difficult to compute. This implementation uses native Python arrays.

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
        X_real = array('f', [ 0. ] * N)        # Accumulate the results;
        X_imag = array('f', [ 0. ] * N)
        W_real = math.cos(-2*math.pi/N)        # Twiddle factors;
        W_imag = math.sin(-2*math.pi/N)
        Wj_real = 1.
        Wj_imag = 0.
        for j in range(N1):                    # Compute every subsequence of size N2;

            Xj = recursive_nfft(x[j::N1])
            Xj_real = array('f', [ xr.real for xr in Xj ])
            Xj_imag = array('f', [ xi.imag for xi in Xj ])
            Wkj_real = 1.
            Wkj_imag = 0.
            for k in range(N):

                # Recombine results:
                aux_real = Xj_real[k%N2]*Wkj_real - Xj_imag[k%N2]*Wkj_imag
                aux_imag = Xj_real[k%N2]*Wkj_imag + Xj_imag[k%N2]*Wkj_real
                X_real[k] += aux_real
                X_imag[k] += aux_imag

                # Update twiddle factors:
                aux_real = Wkj_real*Wj_real - Wkj_imag*Wj_imag
                aux_imag = Wkj_real*Wj_imag + Wkj_imag*Wj_real
                Wkj_real = aux_real
                Wkj_imag = aux_imag

            aux_real = Wj_real*W_real - Wj_imag*W_imag
            aux_imag = Wj_real*W_imag + Wj_imag*W_real
            Wj_real = aux_real
            Wj_imag = aux_imag

        return [ Xr + 1j*Xi for Xr, Xi in zip(X_real, X_imag) ]


####################################################################################################
# Main program:
if __name__ == "__main__":
    pass
