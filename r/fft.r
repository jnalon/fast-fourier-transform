####################################################################################################
# Fast Fourier Transform -- R Version
# This version implements Cooley-Tukey algorithm for powers of 2 only.
#
# José Alexandre Nalon
####################################################################################################
# Since R is an interpreted language, all you have to do is to invoque the interpreter to run
# this program:
#
# $ Rscript fft.r


####################################################################################################
# Definitions:
REPEAT = 50                                    # Number of executions to compute average time;


####################################################################################################
# Measure execution time through repeated calls to a (Fast) Fourier Transform function.
#
# Parameters:
#   f
#     Function to be called;
#   size
#     Power of two of the size of the vector on which the transform will be applied;
#   repeats
#     Number of times the function will be called. Defaults to REPEAT.
#
# Returns:
#   The average execution time for that function with a vector of the given size.
####################################################################################################
time_it <- function(f, size, repeats=REPEAT)
{
    x <- seq(0, size) + 0i                     # Generate a vector;
    t0 <- Sys.time()                           # Start a timer;
    for (j in 0:repeats) {                     # Repeated calls;
        X <- f(x)
    }
    return ((Sys.time() - t0) / repeats)       # Compute average;
}


####################################################################################################
# Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
#
# Parameters:
#   x
#     The vector of which the DFT will be computed. Given the nature of the implementation, there is
#     no restriction on the size of the vector, although it will almost always be called with a
#      power of two size to give a fair comparison.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
direct_ft <- function(x)
 {
    N <- length(x)                             # Length of the vector;
    X <- rep(0+0i, N)                          # Accumulate the results;
    W <- exp(-2i*pi/N)                         # Twiddle factors;
    Wk <- 1.0
    for (k in 1:N) {                           # Compute the kth coefficient;
        Wkn <- 1.0
        for (n in 1:N) {                       #   Operate the summation;
            X[k] <- X[k] + x[n]*Wkn            #     Compute every term;
            Wkn = Wkn * Wk                     # Update twiddle factors;
        }
        Wk = Wk * W
    }
    return (X)
}


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
recursive_fft <- function(x)
{
    if (length(x) == 1) {                                      # A length-1 vector is its own FT;
        return (x)
    } else {
        N <- length(x)                                         # Length of the vector;
        xe <- x[seq(1, N, by=2)]                               # Split into even and odd samples;
        xo <- x[seq(2, N, by=2)]
        Xe <- recursive_fft(xe)                                # Transform of even samples;
        Xo <- recursive_fft(xo)                                # Transform of odd samples;
        W <- exp(-2i*pi/N) ** seq(0, as.integer(N/2)-1)        # Twiddle factors;
        WXo <- W * Xo                                          # Repeated computation;
        X <- c(Xe + WXo, Xe - WXo)                             # Recombine results;
        return (X)
    }
}



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
bit_reverse <- function(k, r)
{
    l <- 0                                     # Accumulate the results;
    for (i in 1:r) {                           # Loop on every bit;
        l <- 2*l + (k %% 2)
        k = k %/% 2                            # Test next bit;
    }
    return (l)
}


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
iterative_fft <- function(x)
{
    N <- length(x)                             # Length of vector;
    r <- as.integer(log2(N))                   # Number of bits;
    X <- x + 0i                                # Accumulate the results;
    for (k in seq(0, N-1)) {
        l <- bit_reverse(k, r)                 # Reorder the vector according to the
        X[l+1] <- x[k+1]                       #   bit-reversed order;
    }

    step <- 1                                  # Auxililary for computation of twiddle factors;
    for (k in 1:r) {
        for (l in seq(0, N-1, by=2*step)) {
            W <- exp(-1i*pi/step)              # Twiddle factors;
            Wkn <- 1.0
            for (n in seq(0, step-1)) {
                p <- l + n + 1
                q <- p + step
                X[q] <- X[p] - Wkn*X[q]        # Recombine results;
                X[p] <- 2*X[p] - X[q]
                Wkn <- Wkn * W                 # Update twiddle factors;
            }
        }
        step <- 2*step
    }
    return (X)
}



####################################################################################################
# Main program:

# Try it with vectors with size ranging from 32 to 1024 samples:
print("+---------+---------+---------+---------+---------+---------+---------+")
print("|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  | Intern. |")
print("+---------+---------+---------+---------+---------+---------+---------+")

# Compute the average execution time:
for (r in 5:10) {

    # Compute the average execution time:
    n = 2**r
    dtime <- time_it(direct_ft, n, REPEAT)
    rtime <- time_it(recursive_fft, n, REPEAT)
    itime <- time_it(iterative_fft, n, REPEAT)
    intime <- time_it(fft, n, REPEAT)

    # Print the results:
    print(sprintf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |",
                  n, n**2, r*n, dtime, rtime, itime, intime))
}

print("+---------+---------+---------+---------+---------+---------+---------+")
