####################################################################################################
# Fast Fourier Transform -- Julia Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# Since Julia is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ julia fft.jl


####################################################################################################
# Import needed modules:
import Printf                                  # Print messages;
import FFTW                                    # Internal Fast Fourier Transform, for comparison;


####################################################################################################
# Definitions:
const REPEAT = 500                             # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function:
"""
    time_it(f, size, repeat=REPEAT)

Average execution time of a function.

This function calls a Fast Fourier Transform function `f` on a vector of `size` elements a certain
number of times (given by `repeat`), and returns the average execution time of the function.
"""
function time_it(f, size, repeat=REPEAT)
    x = Array{Float64}(0:size-1)               # Generate a vector;
    for j = 1:REPEAT                           # Repeated calls;
        f(x)
    end
end


####################################################################################################
"""
    direct_ft(x)

Fourier Transform of a vector `x` by direct definition.

Compute the Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
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

Compute the Fast Fourier Transform using a recursive decimation in time algorithm. This has
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
"""
    bit_reverse(r, k)

Bit-reversed version of an integer number `k` with respect to `r` bits.
"""
function bitreverse(k, r)
    l = 0                                      # Accumulate the results;
    for i = 1:r                                # Loop on every bit;
        l = (l << 1) + (k & 1)                 # Test less signficant bit and add;
        k = k >> 1                             # Test next bit;
    end
    return l
end

"""
    iterative_fft(x)

Fast Fourier Transform of a vector `x` with an iterative decimation-in-time algorithm.

Compute the Fast Fourier Transform using an iterative in-place decimation in time algorithm. This
has O(N log_2(N)) complexity, and since there are less function calls, it will probably be
marginally faster than the recursive versions. The function returns a complex-number vector of the
same size, with the coefficients of the DFT.
"""
function iterative_fft(x)

    N = length(x)                              # Length of vector;
    r = round(Int, log(N)/log(2))              # Number of bits;
    X = complex(x)                             # Accumulate the results;
    for k = 0:N-1
        l = bitreverse(k, r)                   # Reorder the vector according to the
        X[l+1] = x[k+1]                        #   bit-reversed order;
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

    # Start by printing the table with time comparisons:
    println("+---------+---------+---------+---------+---------+---------+---------+")
    println("|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  | Interna |")
    println("+---------+---------+---------+---------+---------+---------+---------+")

    # Try it with vectors with size ranging from 32 to 1024 samples:
    for r = 5:10

        # Compute the average execution time:
        n = 2^r
        dtime = @elapsed time_it(direct_ft, n, REPEAT)
        rtime = @elapsed time_it(recursive_fft, n, REPEAT)
        itime = @elapsed time_it(iterative_fft, n, REPEAT)
        intime = @elapsed time_it(FFTW.fft, n, REPEAT)
        n = 2^r

        # Print the results:
        Printf.@printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n",
                        n, n^2, r*n, dtime, rtime, itime, intime)
    end

    println("+---------+---------+---------+---------+---------+---------+---------+")
end

main()
