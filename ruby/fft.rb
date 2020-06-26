####################################################################################################
# Fast Fourier Transform -- Ruby Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
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
# Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
# complexity.
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
    if x.length == 1                           # A length-1 vector is its own FT;
        return x
    else
        nx = x.length                          # Length of the vector;
        n2 = nx / 2

        xe = [ Complex(0.0, 0.0) ] * n2        # Allocate memory for computation;
        xo = [ Complex(0.0, 0.0) ] * n2
        tX = [ Complex(0, 0) ] * nx

        n2.times { |k|                         # Split even and odd samples;
            xe[k] = x[2*k]
            xo[k] = x[2*k+1]
        }
        tXe = recursive_fft(xe)                # Transform of even samples;
        tXo = recursive_fft(xo)                # Transform of odd samples;

        w = CMath.exp((-2.0i*Math::PI/nx))     # Twiddle factors;
        wk = Complex(1.0, 0.0)
        n2.times { |k|
            tX[k] = tXe[k] + wk*tXo[k]
            tX[k+n2] = tXe[k] - wk*tXo[k]      # Recombine results;
            wk = wk * w                        # Update twiddle factors;
        }
        return tX
    end
end


####################################################################################################
# Bit-reversed version of an integer number.
#
# Parameters:
#   k
#     The number to be bit-reversed;
#   r
#     The number of bits to take into consideration when reversing.
#
# Returns:
#   The number k, bit-reversed according to integers with r bits.
####################################################################################################
def bitreverse(k, r)
    l = 0                                      # Accumulate the results;
    r.times { |i|                              # Loop on every bit;
        l = (l << 1) + (k & 1)                 # Test less signficant bit and add;
        k = (k >> 1)                           # Test next bit;
    }
    return l
end


####################################################################################################
# Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
# O(N log_2(N)) complexity, and since there are less function calls, it will probably be marginally
# faster than the recursive versions. It uses native Python lists.
#
# Parameters:
#   x
#     The vector of which the FFT will be computed. This should always be called with a vector of a
#     power of two length, or it will fail. No checks on this are made.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
def iterative_fft(x)
    nx = x.length                                      # Length of vector;
    r = (Math.log(nx) / Math.log(2)).to_i              # Number of bits;
    tX = [ Complex(0, 0) ] * nx                        # Accumulate the results;
    nx.times { |k|
        l = bitreverse(k, r)                           # Reorder the vector according to the
        tX[l] = x[k]                                   #   bit-reversed order;
    }

    step = 1                                           # Auxililary for computation of twiddle factors;
    r.times { |k|
        l = 0
        while l < nx do
            w = CMath.exp(-1.0i*Math::PI/step)         # Twiddle factors;
            wkn = 1.0
            step.times { |n|
                p = l + n
                q = p + step
                tX[q] = tX[p] - wkn*tX[q]              # Recombine results;
                tX[p] = 2*tX[p] - tX[q]
                wkn = wkn * w                          # Update twiddle factors;
            }
            l = l + 2*step
        end
        step = 2 * step
    }
    return tX
end


####################################################################################################
# Main function:
def main()

    # Try it with vectors with size ranging from 32 to 1024 samples:
    puts "+---------"*6 + "+"
    puts "|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |"
    puts "+---------"*6 + "+"

    # Compute the average execution time:
    5.upto(10) { |r|

        # Compute the average execution time:
        n = 2**r
        dtime = time_it(method(:direct_ft), n, REPEAT)
        rtime = time_it(method(:recursive_fft), n, REPEAT)
        itime = time_it(method(:iterative_fft), n, REPEAT)

        # Print the results:
        puts '| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |' % [ n, n**2, r*n, dtime, rtime, itime ]
    }
    puts "+---------"*6 + "+"
end

main()
