####################################################################################################
# Fast Fourier Transform -- Julia Version
# This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
#
# José Alexandre Nalon
####################################################################################################
# Since Julia is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ julia anyfft.jl


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
        for n = 1:N                            #   Operate the summation;
            X[k] = X[k] + x[n]*Wkn             #     Compute every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        end
        Wk = Wk * W
    end
    return X
end


####################################################################################################
"""
    factor(n)

The smallest prime factor of a given number `n`. If `n` is prime itself, this value is returned.
"""
function factor(n)
    rn = ceil(Int, sqrt(n))                    # Search up to the square root of the number;
    for i = 2:rn+1
        if n%i == 0                            # When remainder is zero, factor is found;
            return i
        end
    end
    return n
end


"""
    recursive_fft(x)

Fast Fourier Transform of a vector `x` with a recursive decimation-in-time algorithm.

Compute the Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
complexity than the direct FT, though the exact value is difficult to compute. The function returns
a complex-number vector of the same size, with the coefficients of the DFT.
"""
function recursive_fft(x)
    N = length(x)
    N1 = factor(N)                             # Find the smallest factor of the vector length;
    if N1 == N                                 # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else
        N2 = floor(Int, N / N1)                # Decompose in two factors, N1 being prime;
        X = complex(zeros(N))                  # Accumulate the results;
        W = exp(-2im*pi/N)                     # Twiddle factors;
        Wj = 1.0
        for j = 1:N1                           # Compute every subsequence of size N2;
            Xj = recursive_fft(x[j:N1:end])
            Wkj = 1.0
            for k = 1:N
                k2 = (k-1) % N2 + 1
                X[k] = X[k] + Xj[k2] * Wkj     # Recombine results;
                Wkj = Wkj * Wj                 # Update twiddle factors;
            end
            Wj = Wj * W
        end
        return X
    end
end


####################################################################################################
# Programa Principal
function main()

    # Start printing the table with time comparisons:
    println("+---------+---------+---------+---------+---------+")
    println("|    N    |   N^2   | Direct  | Recurs. | Intern. |")
    println("+---------+---------+---------+---------+---------+")

    # Try it with vectors with the given sizes:
    sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]
    for size = sizes

        # Compute the average execution time:
        dtime = @elapsed time_it(direct_ft, size, REPEAT)
        rtime = @elapsed time_it(recursive_fft, size, REPEAT)
        intime = @elapsed time_it(FFTW.fft, size, REPEAT)

        # Print the results:
        Printf.@printf("| %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                        size, size^2, dtime, rtime, intime)
    end

    println("+---------+---------+---------+---------+---------+")
end

main()
