# -*- coding: utf-8 -*-
####################################################################################################
# Fast Fourier Transform -- Python 3 Version using NumPy arrays
# This version implements Cooley-Tukey algorithm for powers of 2 and composite numbers
####################################################################################################
# This is a function only library, and executing this file won't give any results. This uses NumPy
# arrays, so you probably won't get any gains by running it with other Python implementatios.


####################################################################################################
# Import needed modules:
import numpy as np                             # Deal with arrays;


####################################################################################################
# Direct FT:
def direct_ft(x):
    """
    Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
    complexity. This implementation uses NumPy arrays.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    x = np.array(x)
    N = len(x)                                 # Length of the vector;
    X = np.zeros(x.shape, dtype=complex)       # Accumulate the results;
    W = np.exp(-2j*np.pi/N)                    # Twiddle factors;
    Wk = 1.
    for k in range(0, N):                      # Compute the kth coefficient;
        Wkn = 1.
        for n in range(0, N):                  #   Operate the summation;
            X[k] = X[k] + x[n]*Wkn             #     Compute every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        Wk = Wk * W
    return X


def simplified_ft(x):
    """
    Discrete Fourier Transform directly from the definition. This implementation uses the power of
    broadcasting from NumPy to make the equations even simpler.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    x = np.array(x)
    N = len(x)
    X = np.zeros(x.shape, dtype=complex)
    W = np.exp(-2j*np.pi/N)
    for k in range(0, N):
        n = np.arange(N)
        Wkn = W ** (k*n)
        X[k] = np.sum(x * Wkn)
    return X


def matrix_ft(x):
    """
    This implementation uses NumPy capabilities to the maximum, by creating a matrix with Fourier
    Transform coefficients (the transform kernel) and matrix multiplying it to the input sequence.
    This is very fast thanks to the internal implementation of matrix product in NumPy, but it's
    still a O(N^2) complexity implementation.

    :Parameters:
      x
        The vector of which the DFT will be computed. Given the nature of the implementation, there
        is no restriction on the size of the vector, although it will almost always be called with
        a power of two size to give a fair comparison.

    :Returns:
      A complex-number vector of the same size, with the coefficients of the DFT.
    """
    x = np.array(x)
    N = len(x)
    k, n = np.meshgrid(np.arange(N), np.arange(N))
    W = np.exp(-2j*k*n*np.pi/N)
    return np.dot(W, x)


####################################################################################################
# Recursive Decimation-in-time FFT using NumPy arrays:
def recursive_fft(x):
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
        return np.array(x)
    else:
        x = np.array(x)
        N = len(x)                                     # Length of the vector;
        Xe = recursive_fft(x[0::2])                    # Transform of even samples;
        Xo = recursive_fft(x[1::2])                    # Transform of odd samples;
        W = np.exp(-2j*np.pi/N) ** range(0, N//2)      # Twiddle factors;
        WXo = W * Xo                                   # Repeated computation;
        X = np.hstack((Xe + WXo, Xe - WXo))            # Recombine results;
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
# Iterative Decimation-in-time FFT, using NumPy arrays:
def iterative_fft(x):
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
    x = np.array(x)
    N = len(x)                                 # Length of vector;
    r = int(np.log2(N))                        # Number of bits;
    X = np.array(x, dtype=complex)             # Accumulate the results;
    for k in range(0, N):
        l = bit_reverse(k, r)                  # Reorder the vector according to the
        X[l] = x[k]                            #   bit-reversed order;

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k in range(0, r):
        for l in range(0, N, 2*step):
            W = np.exp(-1j * np.pi / step)      # Twiddle factors;
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
    rn = int(np.ceil(np.sqrt(n)))              # Search up to the square root of the number;
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
    value is difficult to compute. This implementation uses NumPy arrays.

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
        x = np.array(x)
        N2 = N // N1                           # Decompose in two factors, N1 being prime;
        X = np.zeros((N, ), dtype=complex)     # Accumulate the results;
        W = np.exp(-2j*np.pi/N)                # Twiddle factors;
        Wj = 1.
        for j in range(N1):                    # Compute every subsequence of size N2;
            Xj = recursive_nfft(x[j::N1])      # Recursively compute the Fourier Transform;
            Wkj = 1.
            for k in range(N):
                X[k] = X[k] + Xj[k%N2] * Wkj   # Recombine results;
                Wkj = Wkj * Wj                 # Update twiddle factors;
            Wj = Wj * W
        return X


####################################################################################################
# Recursive FFT:
def vec_recursive_nfft(x):
    """
    Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
    complexity than the direct FT, though the exact value is difficult to compute. This
    implementation uses NumPy arrays for conciseness. In this implementation, loops are avoided by
    vectorizing the computation of the twiddle factors.

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
        X = np.zeros((N, ), dtype=complex)     # Accumulate the results;
        k = np.arange(N)
        Wk = np.exp(-2j*np.pi*k/N)             # Twiddle factors;
        Wkj = np.ones(Wk.shape)
        for j in range(N1):                    # Compute every subsequence of size N2;
            Xj = vec_recursive_nfft(x[j::N1])  # Recursively compute the Fourier Transform;
            X = X + Xj[k%N2]*Wkj               # Recombine results;
            Wkj = Wkj * Wk                     # Update twiddle factors;
        return X


####################################################################################################
# Main program:
if __name__ == "__main__":
    pass
