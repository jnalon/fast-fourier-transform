'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This version implements Cooley-Tukey algorithm for powers of 2 and composite numbers.
'' -------------------------------------------------------------------------------------------------

'' Include necessary libraries:
#include "complex.bas"


'' Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
''
'' Parameters:
''   X
''     The vector of which the DFT will be computed. Given the nature of the implementation, there
''     is no restriction on the size of the vector, although it will almost always be called with a
''     power of two size to give a fair comparison;
''   TX
''     The vector that will receive the results of the computation. It needs to be declared by the
''     caller. The result will be computed in place.
sub DirectFT(X() as Complex, TX() as Complex)
    dim Nx as Integer = ubound(x)
    dim W as Complex = Complex.exp(-2*PI/(Nx+1))       '' Initialize twiddle factors;
    dim Wk as Complex = Complex(1, 0)
    dim Wkn as Complex
    dim k as Integer, n as Integer
    redim TX(Nx)
    for k = 0 to Nx
        TX(k) = Complex()                              '' Accumulate the results;
        Wkn = Complex(1, 0)                            '' Initialize twiddle factors;
        for n = 0 to Nx
            TX(k) = TX(k) + Wkn * X(n)
            Wkn = Wkn * Wk                             '' Update twiddle factors;
        next n
        Wk = Wk * W
    next k
end sub


'' Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
'' complexity.
''
'' Parameters:
''   X
''     The vector of which the FFT will be computed. This should always be called with a vector of a
''     power of two length, or it will fail. No checks on this are made.
''   TX
''     The vector that will receive the results of the computation. It needs to be declared by the
''     caller. The result will be computed in place.
sub RecursiveFFT(X() as Complex, TX() as Complex)
    dim Nx as Integer = ubound(x)                      '' A length-1 vector is its own FT;
    redim TX(Nx)
    if Nx = 0 then
        TX(0) = X(0)
    else
        dim N2 as Integer = (Nx+1)\2 - 1
        dim Xe(N2) as Complex                          '' Allocate memory for computation;
        dim Xo(N2) as Complex
        dim TXe() as Complex
        dim TXo() as Complex
        dim k as Integer
        for k = 0 to N2                                '' Split even and odd samples;
            Xe(k) = X(2*k)
            Xo(k) = X(2*k + 1)
        next k
        RecursiveFFT(Xe(), TXe())                      '' Transform of even samples;
        RecursiveFFT(Xo(), TXo())                      '' Transform of odd samples;

        dim W as Complex = Complex.exp(-2*PI/(Nx+1))   '' Twiddle factors;
        dim Wk as Complex = Complex(1, 0)
        dim WkXok as Complex
        for k = 0 to N2
            WkXok = Wk * TXo(k)
            TX(k) = TXe(k) + WkXok                     '' Recombine results;
            TX(k+N2+1) = TXe(k) - WkXok
            Wk = Wk * W                                '' Update twiddle factors;
        next k
    end if
end sub


'' Bit-reversed version of an integer number.
''
'' Parameters:
''   k
''     The number to be bit-reversed;
''   r
''     The number of bits to take into consideration when reversing.
''
'' Returns:
''   The number k, bit-reversed according to integers with r bits.
private function BitReverse(k as Integer, r as Integer) as Integer
    dim l as Integer = 0                               '' Accumulate the results;
    dim i as Integer
    for i = 1 to r                                     '' Loop on every bit;
        l = (2 * l) + (k mod 2)                        '' Test less signficant bit and add;
        k = k \ 2                                      '' Test next bit;
    next i
    return l
end function


'' Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
'' O(N log_2(N)) complexity, and since there are less function calls, it will probably be marginally
'' faster than the recursive versions.
''
'' Parameters:
''   X
''     The vector of which the FFT will be computed. This should always be called with a vector of a
''     power of two length, or it will fail. No checks on this are made.
''   TX
''     The vector that will receive the results of the computation. It needs to be declared by the
''     caller. The result will be computed in place.
sub IterativeFFT(X() as Complex, TX() as Complex)
    dim Nx as Integer = ubound(X)
    dim r as Integer = int(log(Nx+1)/log(2))           '' Number of bits;
    dim W as Complex, Wkn as Complex
    dim k as Integer, l as Integer, n as Integer, p as Integer, q as Integer, stp as Integer

    redim TX(Nx)
    for k = 0 to Nx
        l = BitReverse(k, r)                           '' Reorder the vector according to the
        TX(l) = X(k)                                   ''   bit-reversed order;
    next k

    stp = 1                                            '' Computation of twiddle factors;
    for k = 1 to r
        for l = 0 to Nx-1 step 2*stp
            W = Complex.exp(-PI/stp)                   '' Twiddle factors;
            Wkn = Complex(1, 0)
            for n = 0 to stp-1
                p = l + n
                q = p + stp
                TX(q) = TX(p) - Wkn * TX(q)            '' Recombine results;
                TX(p) = 2 * TX(p) - TX(q)
                Wkn = Wkn * W                          '' Update twiddle factors;
            next n
        next l
        stp = 2 * stp
    next k
end sub


'' Smallest prime factor of a given number. If the argument is prime itself, then it is the return
'' value.
''
'' Parameters:
''   n
''     Number to be inspected.
''
'' Returns:
''   The smallest prime factor, or the number itself if it is already a prime.
private function Factor(N as Integer) as Integer
    dim rn as Integer = int(sqr(N)) + 1
    dim i as Integer
    for i = 2 to rn+1
        if N mod i = 0 then
            return i
        end if
    next i
    return N
end function


sub RecursiveNFFT(X() as Complex, TX() as Complex)
    dim Nx as Integer = ubound(X)
    dim N1 as Integer, j as Integer, n as Integer, k as Integer

    N1 = Factor(Nx+1)
    if N1 = Nx + 1 then
        DirectFT(X(), TX())
    else
        dim Xj() as Complex, TXj() as Complex
        dim W as Complex, Wj as Complex, Wkj as Complex
        dim N2 as Integer = (Nx + 1) \ N1
        redim Xj(N2-1)
        redim TXj(N2-1)
        W = Complex.exp(-2*PI/(Nx+1))
        Wj = Complex(1, 0)
        for j = 0 to N1 - 1
            for n = 0 to N2 - 1
                Xj(n) = X(n*N1+j)
                TXj(n) = Complex()
            next n
            RecursiveNFFT(Xj(), TXj())
            Wkj = Complex(1, 0)
            for k = 0 to Nx
                TX(k) = TX(k) + Wkj * TXj(k mod N2)
                Wkj = Wkj * Wj
            next k
            Wj = Wj * W
        next j
    end if
end sub


