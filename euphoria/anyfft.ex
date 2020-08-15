----------------------------------------------------------------------------------------------------
-- Fast Fourier Transform -- Euphoria 4 Version
-- This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
--
-- José Alexandre Nalon
----------------------------------------------------------------------------------------------------
-- Euphoria can be interpreted or compiled. Since this program is only a simple implementation to
-- test the basic capabilities of the language, it doesn't need much to run. Just type in the
-- command line:
--
-- $ eui anyfft.ex
----------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------
-- Import needed modules:
include std/math.e                             -- Math functions;


----------------------------------------------------------------------------------------------------
-- Definitions:
constant REPEAT = 50                           -- Number of executions to compute average time;


----------------------------------------------------------------------------------------------------
-- Small library to deal with complex numbers:
enum REAL, IMAG
constant PI = 3.1415926535897931

-- The Complex data structure:
type complex(sequence z)
    return length(z) = 2 and atom(z[1]) and atom(z[2])
end type

-- Product of complex numbers:
function cmul(complex z1, complex z2)
    return { z1[REAL]*z2[REAL] - z1[IMAG]*z2[IMAG],
             z1[REAL]*z2[IMAG] + z1[IMAG]*z2[REAL] }
end function

-- Convenience function to compute exponentials:
function cexp(complex z)
    atom r = power(2.7182818284590451, z[REAL])
    atom theta = z[IMAG]
    return { r*cos(theta), r*sin(theta) }
end function


----------------------------------------------------------------------------------------------------
-- Auxiliary procedure: complex_show
--   Pretty printing of an array of complex numbers, used to inspect results.
--
-- Parameters:
--   x
--     A vector of complex numbers, according to the definition above;
----------------------------------------------------------------------------------------------------
procedure complex_show(sequence x)
    for j = 1 to length(x) do
        printf(1, "(%7.4f, %7.4f)\n", { x[j][1], x[j][2] })
    end for
end procedure


----------------------------------------------------------------------------------------------------
-- Measure execution time through repeated calls to a (Fast) Fourier Transform function.
--
-- Parameters:
--   f
--     Function to be called;
--   size
--     Power of two of the size of the vector on which the transform will be applied;
--   repeats
--      Number of times the function will be called. Defaults to REPEAT.
--
-- Returns:
--   The average execution time for that function with a vector of the given size.
----------------------------------------------------------------------------------------------------
function time_it(atom f, integer size, integer repeats=REPEAT)
    sequence x = repeat({ 0.0, 0.0 }, size)    -- Generate a complex vector;
    for j = 1 to size do
        x[j] = { j-1, 0.0 }
    end for

    atom t0 = time()                           -- Start a timer;
    for j = 1 to repeats do                    -- Repeated calls;
        sequence X = call_func(f, { x })
    end for
    return (time() - t0) / repeats             -- End timer and compute average times;
end function


----------------------------------------------------------------------------------------------------
-- Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
--
-- Parameters:
--   x
--     The vector of which the DFT will be computed. Given the nature of the implementation, there
--     is no restriction on the size of the vector, although it will almost always be called with a
--     a power of two size to give a fair comparison.
--
-- Returns:
--   A complex-number vector of the same size, with the coefficients of the DFT.
----------------------------------------------------------------------------------------------------
function direct_ft(sequence x)
    integer N = length(x)                      -- Length of the vector;
    complex W = cexp({ 0.0, -2*PI/N })         -- Twiddle factors;

    sequence X = repeat({ 0.0, 0.0 }, N)       -- Initialize the results;
    complex Wk = { 1.0, 0.0 }                  -- Initialize twiddle factor;
    for k = 1 to N do                          -- Compute kth coefficient;
        complex Wkn = { 1.0, 0.0 }             --   Initialize twiddle factors;
        for n = 1 to N do                      --   Operate the summation;
            X[k] = X[k] + cmul(x[n], Wkn)      --     Compute every term;
            Wkn = cmul(Wkn, Wk)                --     Update twiddle factor;
        end for
        Wk = cmul(Wk, W)                       -- Update twiddle factor;
    end for
    return X
end function


----------------------------------------------------------------------------------------------------
-- Smallest prime factor of a given number. If the argument is prime itself, then it is the return
-- value.
--
-- Parameters:
--   n
--     Number to be inspected.
--
-- Returns:
--   The smallest prime factor, or the number itself if it is already a prime.
----------------------------------------------------------------------------------------------------
function factor(integer n)
    integer rn = floor(n/2)                    -- Search up to half the number;
    for i=2 to rn do
        if mod(n, i) = 0 then                  -- If remainder is zero, a factor is found;
            return i
        end if
    end for
    return n
end function


----------------------------------------------------------------------------------------------------
-- Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
-- complexity than the direct FT, though the exact value is difficult to compute.
--
-- Parameters:
--   x
--     The vector of which the FFT will be computed. Its length must be a composite number, or else
--     the computation will be defered to the direct FT, and there will be no efficiency gain.
--
-- Returns:
--   A complex-number vector of the same size, with the coefficients of the DFT.
----------------------------------------------------------------------------------------------------
function recursive_fft(sequence x)
    integer N = length(x)
    integer N1 = factor(N)                             -- Smallest prime factor of length;

    if N1 = N then                                     -- If the length is prime itself,
        return direct_ft(x)
    else
        integer N2 = floor(N/N1)                       -- Decompose in two factors, N1 being prime;

        sequence X = repeat({ 0.0, 0.0 }, N)           -- Allocate memory for computation
        sequence xj = repeat({ 0.0, 0.0 }, N2)         -- Allocate memory for subsequences

        complex W = cexp({ 0.0, -2.0*PI/N })           -- Twiddle factor;
        complex Wj = { 1.0, 0.0 }
        for j = 1 to N1 do                             -- Compute every subsequence of size N2;
            for n = 0 to N2-1 do
                xj[n+1] = x[n*N1+j]
            end for
            sequence Xj = recursive_fft(xj)
            complex Wkj = { 1.0, 0.0 }
            for k = 1 to N do
                X[k] = X[k] + cmul(Xj[mod(k-1, N2)+1], Wkj)    -- Recombine results;
                Wkj = cmul(Wkj, Wj)                            -- Update twiddle factors;
            end for
            Wj = cmul(Wj, W)
        end for
        return X
    end if
end function


----------------------------------------------------------------------------------------------------
-- Main function:
procedure main()

    sequence SIZES = { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 }

    -- Start by printing the table with time comparisons:
    printf(1, "%s\n", { "+---------+---------+---------+---------+" })
    printf(1, "%s\n", { "|    N    |   N^2   | Direct  | Recurs. |" })
    printf(1, "%s\n", { "+---------+---------+---------+---------+" })

    -- Try it with vectors with the given sizes:
    for i = 1 to length(SIZES) do

        -- Compute the average execution time:
        integer n = SIZES[i]
        atom dtime = time_it(routine_id("direct_ft"), n, REPEAT)
        atom rtime = time_it(routine_id("recursive_fft"), n, REPEAT)

        -- Print the results:
        sequence tup = { n, power(n, 2), dtime, rtime }
        printf(1, "| %7d | %7d | %7.4f | %7.4f |\n", tup)

    end for
    printf(1, "%s\n", { "+---------+---------+---------+---------+" })

end procedure


main()
