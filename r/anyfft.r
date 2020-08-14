####################################################################################################
# Fast Fourier Transform -- R Version
# This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
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
    (Sys.time() - t0) / repeats                # Compute average;
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
    X
}


####################################################################################################
# Smallest prime factor of a given number. If the argument is prime itself, then it is the return
# value.
#
# Parameters:
#   n
#     Number to be inspected.
#
# Returns:
#   The smallest prime factor, or the number itself if it is already a prime.
####################################################################################################
factor <- function(n)
{
    rn <- ceiling(sqrt(n))                     # Search up to the square root of the number;
    for (i in seq(2, rn)) {
        if (n%%i == 0) {                       # When remainder is zero, factor is found;
            return (i)
        }
    }
    return (n)
}


####################################################################################################
# Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller complexity
# than the direct FT, though the exact value is difficult to compute.
#
# Parameters:
#   x
#     The vector of which the FFT will be computed. Its length must be a composite number, or else
#     the computation will be defered to the direct FT, and there will be no efficiency gain.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
recursive_fft <- function(x)
{
    N <- length(x)                             # Length of the vector;
    N1 <- factor(N)                            # Find the smallest factor of the vector length;
    if (N1 == N) {                             # If the length is prime itself,
        return (direct_ft(x))                  #    the transform is given by the direct form;
    } else {
        N2 <- N %/% N1                         # Decompose in two factors, N1 being prime;
        X <- rep(0+0i, N)                      # Accumulate the results;
        W <- exp(-2i*pi/N)                     # Twiddle factors;
        Wj <- 1.0
        for (j in 1:N1) {                      # Compute every subsequence of size N2;
            xj = x[seq(j, N, by=N1)]
            Xj = recursive_fft(xj)
            Wkj <- 1.0
            for (k in 1:N) {
                k2 <- (k-1) %% N2              # Recombine results;
                X[k] <- X[k] + Xj[k2+1] * Wkj  # Update twiddle factors;
                Wkj <- Wkj * Wj
            }
            Wj <- Wj * W
        }
        return (X)
    }
}


####################################################################################################
# Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller complexity
# than the direct FT, though the exact value is difficult to compute. In this implementation, loops
# are avoided by vectorizing the computation of the twiddle factors.
#
# Parameters:
#   x
#     The vector of which the FFT will be computed. Its length must be a composite number, or else
#     the computation will be defered to the direct FT, and there will be no efficiency gain.
#
# Returns:
#   A complex-number vector of the same size, with the coefficients of the DFT.
####################################################################################################
vec_recursive_fft <- function(x)
{
    N <- length(x)                             # Length of the vector;
    N1 <- factor(N)                            # Find the smallest factor of the vector length;
    if (N1 == N) {                             # If the length is prime itself,
        return (direct_ft(x))                  #    the transform is given by the direct form;
    } else {
        N2 <- N %/% N1                         # Decompose in two factors, N1 being prime;
        X <- rep(0+0i, N)                      # Accumulate the results;
        k <- seq(0, N-1)
        Wk <- exp(-2i*pi*k/N)                  # Twiddle factors;
        Wkj <- rep(1+0i, N)
        for (j in 1:N1) {                      # Compute every subsequence of size N2;
            xj = x[seq(j, N, by=N1)]           # Split subsequences;
            Xj <- vec_recursive_fft(xj)        # Recursively compute the Fourier Transform;
            k2 <- k %% N2
            X = X + Xj[k2+1] * Wkj             # Recombine results;
            Wkj = Wkj * Wk                     # Update twiddle factors;
        }
        return (X)
    }
}



####################################################################################################
# Main program:

SIZES = c( 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 )

# Try it with vectors with size ranging from 32 to 1024 samples:
print("+---------+---------+---------+---------+---------+---------+")
print("|    N    |   N^2   | Direct  | Recurs. | VecRec. | Intern. |")
print("+---------+---------+---------+---------+---------+---------+")

# Compute the average execution time:
for (n in SIZES) {

    # Compute the average execution time:
    dtime <- time_it(direct_ft, n, REPEAT)
    rtime <- time_it(recursive_fft, n, REPEAT)
    vtime <- time_it(vec_recursive_fft, n, REPEAT)
    intime <- time_it(fft, n, REPEAT)

    # Print the results:
    print(sprintf("| %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |",
                  n, n**2, dtime, rtime, vtime, intime))
}

print("+---------+---------+---------+---------+---------+---------+")
