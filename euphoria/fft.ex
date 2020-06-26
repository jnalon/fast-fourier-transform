----------------------------------------------------------------------------------------------------
-- Fast Fourier Transform -- Euphoria 4 Version
-- This version implements Cooley-Tukey algorithm for powers of 2 only.
--
-- José Alexandre Nalon
----------------------------------------------------------------------------------------------------
-- Euphoria can be interpreted or compiled. Since this program is only a simple implementation to
-- test the basic capabilities of the language, it doesn't need much to run. Just type in the
-- command line:
--
-- $ eui fft.ex
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
-- Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
-- complexity.
--
-- Parameters:
--   x
--     The vector of which the FFT will be computed. This should always be called with a vector of a
--     power of two length, or it will fail. No checks on this are made.
--
-- Returns:
--   A complex-number vector of the same size, with the coefficients of the DFT.
----------------------------------------------------------------------------------------------------
function recursive_fft(sequence x)
    integer N = length(x)

    if length(x) = 1 then                              -- A length-1 vector is its own FT;
        return x
    else
        integer N2 = floor(N/2)

        sequence X = repeat({ 0.0, 0.0 }, N)           -- Allocate memory for computation;
        sequence xe = repeat({ 0.0, 0.0 }, N2)
        sequence xo = repeat({ 0.0, 0.0 }, N2)

        for k = 1 to N2 do                             -- Split even and odd samples;
            xe[k] = x[2*k-1]
            xo[k] = x[2*k]
        end for
        sequence Xe = recursive_fft(xe)                -- Transform of even samples;
        sequence Xo = recursive_fft(xo)                -- Transform of odd samples;

        complex W = cexp({ 0.0, -2.0*PI/N })           -- Twidddle factors;
        complex Wk = { 1.0, 0.0 }
        complex wkn_Xo
        for k = 1 to N2 do
            wkn_Xo = cmul(Wk, Xo[k])                   -- Recombine results;
            X[k] = Xe[k] + wkn_Xo
            X[k+N2] = Xe[k] - wkn_Xo
            Wk = cmul(Wk, W)                           -- Update twiddle factors;
        end for
        return X
    end if
end function


----------------------------------------------------------------------------------------------------
-- Bit-reversed version of an integer number.
--
-- Parameters:
--   k
--     The number to be bit-reversed;
--   r
--     The number of bits to take into consideration when reversing.
--
-- Returns:
--   The number k, bit-reversed according to integers with r bits.
----------------------------------------------------------------------------------------------------
function bit_reverse(integer k, integer r)
    integer l = 0                              -- Accumulate the results;
    for i = 1 to r do                          -- Loop on every bit;
        l = 2*l + mod(k, 2)                    -- Test less signficant bit and add;
        k = floor(k/2)                         -- Test next bit;
    end for
    return l
end function


----------------------------------------------------------------------------------------------------
-- Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
-- O(N log_2(N)) complexity, and since there are less function calls, it will probably be marginally
-- faster than the recursive versions.
--
-- Parameters:
--   x
--     The vector of which the FFT will be computed. This should always be called with a vector of a
--     power of two length, or it will fail. No checks on this are made.
--
-- Returns:
--   A complex-number vector of the same size, with the coefficients of the DFT.
----------------------------------------------------------------------------------------------------
function iterative_fft(sequence x)
    integer N = length(x)                              -- Length of the vector;
    integer r = floor(log(N)/log(2))                   -- Number of bits;
    sequence X = x
    for k = 0 to N-1 do
        integer l = bit_reverse(k, r)                  -- Reorder the vector according to the
        X[l+1] = x[k+1]                                --    bit-reversed order;
    end for

    integer step = 1                                   -- Auxiliary for computation of twiddle factors;
    for k = 1 to r do
        for l = 0 to N-1 by 2*step do
            complex W = cexp({ 0.0, -PI/step })        -- Twiddle factors;
            complex Wkn = { 1.0, 0.0 }
            for n = 1 to step do
                integer p = l + n
                integer q = p + step
                X[q] = X[p] - cmul(Wkn, X[q])          -- Recombine results;
                X[p] = X[p] + X[p] - X[q]
                Wkn = cmul(Wkn, W)                     -- Update twiddle factors;
            end for
        end for
        step = 2*step
    end for
    return X
end function


----------------------------------------------------------------------------------------------------
-- Main function:
procedure main()

    -- Start by printing the table with time comparisons:
    printf(1, "%s\n", { "+---------+---------+---------+---------+---------+---------+" })
    printf(1, "%s\n", { "|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |" })
    printf(1, "%s\n", { "+---------+---------+---------+---------+---------+---------+" })

    -- Try it with vectors with size ranging from 32 to 1024 samples:
    for r = 5 to 10 do

        -- Compute the average execution time:
        integer n = power(2, r)
        atom dtime = time_it(routine_id("direct_ft"), n, REPEAT)
        atom rtime = time_it(routine_id("recursive_fft"), n, REPEAT)
        atom itime = time_it(routine_id("iterative_fft"), n, REPEAT)

        -- Print the results:
        sequence tup = { n, power(n, 2), r*n, dtime, rtime, itime }
        printf(1, "| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n", tup)

    end for
    printf(1, "%s\n", { "+---------+---------+---------+---------+---------+---------+" })

end procedure


main()
