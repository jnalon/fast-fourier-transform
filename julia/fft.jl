####################################################################################################
# Fast Fourier Transform -- Julia Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# Since Python is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ julia fft.jl


####################################################################################################
# Import needed modules:
import Printf                                  # Print messages;
import Dates                                   # Time events;
import FFTW                                    # Internal Fast Fourier Transform, for comparison;


####################################################################################################
# Definitions:
const REPEAT = 500                             # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function:
"""
    time_it(f, r, repeat=REPEAT)

Average execution time of a function.

This function calls a Fast Fourier Transform function `f` on a vector of size `2^r` a certain number
of times (given by `repeat`), and returns the average execution time of the function.
"""
function time_it(f, r, repeat=REPEAT)
    x = Array{Float64}(0:2^r-1)                  # Generate a vector;
#     t0 = Dates.now()                           # Starts a timer;
    for j = 1:REPEAT                           # Repeated calls;
        f(x)
    end
#     return float(Dates.now() - t0) / 1000.0 / float(repeat)
end


####################################################################################################
"""
    direct_ft(x)

Fourier Transform of a vector `x` by direct definition.

Computes the Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
complexity. The function returns a complex-number vector of the same size, with the coefficients of
the DFT.
"""
function direct_ft(x)
    N = length(x)                              # Length of the vector;
    X = zeros(Complex, size(x))                # Accumulates the results;
    W = exp(-2im*pi/N)                         # Twiddle factors;
    Wk = 1.0
    for k = 1:N                                # Compute the kth coefficient;
        Wkn = 1
        for n = 1:N                            #   Operates the summation;
            X[k] = X[k] + x[n]*Wkn             #     Computes every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        end
        Wk = Wk * W
    end
    return X
end


####################################################################################################
"""
    recursive_fft(x)

Fast Fourier Transform of a vector `x` with a recursive decimation-in-time algorithm.

Computes the Fast Fourier Transform using a recursive decimation in time algorithm. This has
O(N log_2(N)) complexity. The function returns a complex-number vector of the same size, with the
coefficients of the DFT.
"""
function recursive_fft(x)
    if length(x) == 1                                  # A length-1 vector is its own FT;
        return x
    else
        N = length(x)                                  # Length of the vector;
        N2 = round(Int, N/2)
        Xe = recursive_fft(x[1:2:end])                 # Transform of even samples;
        Xo = recursive_fft(x[2:2:end])                 # Transform of odd samples;
        W = exp.(-2im*pi*Array(0:N2-1)/N)              # Twiddle factors;
        WXo = W .* Xo                                  # Repeated computation;
        X = [ Xe + WXo; Xe-WXo ]                       # Recombine results;
        return X
    end
end


####################################################################################################
function bitreverse(r, k)
    l = 0
    for i = 1:r
        l = (l << 1) + (k & 1)
        k = k >> 1
    end
    return l
end

"""
    iterative_fft(x)

Fast Fourier Transform of a vector `x` with an iterative decimation-in-time algorithm.

Computes the Fast Fourier Transform using an iterative in-place decimation in time algorithm. This
has O(N log_2(N)) complexity, and since there are less function calls, it will probably be
marginally faster than the recursive versions. The function returns a complex-number vector of the
same size, with the coefficients of the DFT.
"""
function iterative_fft(x)

    N = length(x)                              # Length of vector;
    r = round(Int, log(N)/log(2))              # Number of bits;
    X = complex(x)                             # Accumulates the results;
    for k = 0:N-1                              # Reorder the vector according to the
        l = bitreverse(r, k)                   #   bit-reversed order;
        X[l+1] = x[k+1]
    end

    step = 1                                   # Auxililary for computation of twiddle factors;
    for k = 1:r
        for l = 0:2*step:N-1
            W = exp(-1im * pi / step)          # Twiddle factors;
            Wkn = 1.0
            for n = 1:step
                p = l + n
                q = p + step
                X[q] = X[p] - Wkn*X[q]         # Recombine results;
                X[p] = 2*X[p] - X[q]
                Wkn = Wkn * W                  # Update twiddle factors;
            end
        end
        step = 2*step
    end
    return X

end


####################################################################################################
# Programa Principal
function main()
    println("+---------+---------+---------+---------+---------+---------+---------+")
    println("|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  | Interna |")
    println("+---------+---------+---------+---------+---------+---------+---------+")

    for r = 5:10
        dtime = @elapsed time_it(direct_ft, r, REPEAT)
        rtime = @elapsed time_it(recursive_fft, r, REPEAT)
        itime = @elapsed time_it(iterative_fft, r, REPEAT)
        intime = @elapsed time_it(FFTW.fft, r, REPEAT)
        n = 2^r
        Printf.@printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n", n, n^2, r*n, dtime, rtime, itime, intime)
    end

    println("+---------+---------+---------+---------+---------+---------+---------+")
end

main()
