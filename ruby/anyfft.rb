####################################################################################################
# Fast Fourier Transform -- Ruby Version
# This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
#
# José Alexandre Nalon
####################################################################################################
# Since Ruby is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ ruby fft.rb
#
# You might need to install some of the packages imported below.


####################################################################################################
# Import needed modules:
require 'complex'                              # Complex numbers;
require 'cmath'                                # Complex math;
require 'time'                                 # Time measurement;


####################################################################################################
# Definitions:
REPEAT = 50                                    # Number of executions to compute average time;


####################################################################################################
# Auxiliary Function: complex_show
#   Pretty printing of an array of complex numbers, used to inspect results.
#
# Parameters:
#   x
#     A vector of complex numbers, according to the definition above;
####################################################################################################
def complex_show(x)
    x.length.times { |i|
        puts '(%7.4f, %7.4f)' % [ x[i].real, x[i].imag ]
    }
end


####################################################################################################
# Auxiliary Function: time_it
#   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
#
# Parameters:
#   f
#     Function to be called;
#   size
#     Power of two of the size of the vector on which the transform will be applied;
#   repeat
#     Number of times the function will be called. Defaults to REPEAT.
#
# Returns:
#   The average execution time for that function with a vector of the given size.
####################################################################################################
def time_it(f, size, repeat=REPEAT)
    x = Array.new(size) { |e| e = Complex(e, 0) }      # Generate a vector;
    t0 = Time.now                                      # Start a timer;
    repeat.times { y = f.call(x) }                     # Repeated calls;
    return (Time.now - t0) / repeat                    # Compute average;
end


####################################################################################################
# Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
#
# Parameters:
#   x
#     The vector of which the DFT will be computed. Given the nature of the implementation, there is
#     no restriction on the size of the vector, although it will almost always be called with a
#     power of two size to give a fair comparison.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
def direct_ft(x)
    nx = x.length                              # Length of the vector;
    tX = [ Complex(0.0, 0.0) ] * nx            # Accumulate the results;
    w = CMath.exp(-2.0i*Math::PI/nx)           # Twiddle factors;
    wk = 1.0
    nx.times { |k|                             # Compute the kth coefficient;
        wkn = 1.0
        nx.times { |n|                         #   Operate the summation;
            tX[k] = tX[k] + x[n] * wkn         #     Compute every term;
            wkn = wkn * wk                     # Update twiddle factors;
        }
        wk = wk * w
    }
    return tX
end


####################################################################################################
# Auxiliary Function: factor
#   Smallest prime factor of a given number. If the argument is prime itself, then it is the return
#   value.
#
# Parameters:
#   n
#     Number to be inspected.
#
# Returns:
#   The smallest prime factor, or the number itself if it is already a prime.
####################################################################################################
def factor(n)
    rn = Math.sqrt(n).ceil()                   # Search up to the square root of the number;
    2.upto(rn+1) { |i|
        if n%i == 0                            # When remainder is zero, factor is found;
            return i
        end
    }
    return n
end


####################################################################################################
# Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller complexity
# than the direct FT, though the exact value is difficult to compute.
#
# Parameters:
#   x
#     The vector of which the FFT will be computed. This should always be called with a vector of a
#     power of two length, or it will fail. No checks on this are made.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
def recursive_fft(x)
    nx = x.length                              # Length of the vector;
    n1 = factor(nx)                            # Find the smallest factor of the vector length;
    if n1 == nx                                # If the length is prime itself,
        return direct_ft(x)                    #    the transform is given by the direct form;
    else
        n2 = nx / n1                           # Decompose in two factors, N1 being prime;

        xj = [ Complex(0.0, 0.0) ] * n2        # Allocate memory for computation;
        tX = [ Complex(0.0, 0.0) ] * nx        # Accumulate the results;

        w = CMath.exp(-2.0i*Math::PI/nx)       # Twiddle factors;
        wj = 1.0
        n1.times { |j|                         # Compute every subsequence of size N2;
            n2.times { |n|
                xj[n] = x[n*n1+j]
            }
            tXj = recursive_fft(xj)
            wkj = 1.0
            nx.times { |k|
                tX[k] = tX[k] + tXj[k%n2]*wkj  # Recombine results;
                wkj = wkj * wj                 # Update twiddle factors;
            }
            wj = wj * w
        }
        return tX
    end
end


####################################################################################################
# Main function:
def main()

    sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]

    # Try it with vectors with size ranging from 32 to 1024 samples:
    puts "+---------"*4 + "+"
    puts "|    N    |   N^2   | Direct  | Recurs. |"
    puts "+---------"*4 + "+"

    # Compute the average execution time:
    sizes.each { |n|

        # Compute the average execution time:
        dtime = time_it(method(:direct_ft), n, REPEAT)
        rtime = time_it(method(:recursive_fft), n, REPEAT)

        # Print the results:
        puts '| %7d | %7d | %7.4f | %7.4f |' % [ n, n**2, dtime, rtime ]
    }
    puts "+---------"*4 + "+"
end

main()
