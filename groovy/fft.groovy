/**************************************************************************************************
 * Fast Fourier Transform -- Groovy Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 Since Groovy is an interpreted language, all you have to do is to invoque the interpreter to run
 this program:

 $ groovy fft.groovy
 **************************************************************************************************/

/**************************************************************************************************
 Definitions
 **************************************************************************************************/
REPEAT = 50                                    // Number of executions to compute average time;


/**************************************************************************************************
 * Mini-library to deal with complex numbers. While it may be useful to use modules in bigger
 * projects, this one is very small and thorougly tested, and it is simples to just add the
 * code to the main program.
 **************************************************************************************************/
class Complex {

    def re
    def im

    // Constructors:
    Complex(real, imag) {
        re = real
        im = imag
    }

    Complex(c) {
        re = c.re
        im = c.im
    }

    // Add the argument to this, giving the result as a new complex number:
    def plus(c) {
        return new Complex(re + c.re, im + c.im)
    }

    // Subtract the argument from this, giving the result as a new complex number:
    def minus(c) {
        return new Complex(re - c.re, im - c.im)
    }

    // Unitary negation:
    def negative() {
        return new Complex(-re, -im)
    }

    // Multiply this by the argument, giving the result as a new complex number:
    def multiply(c) {
        return new Complex(re*c.re - im*c.im, re*c.im + im*c.re)
    }

    // Complex exponential of an angle:
    static def exp(a) {
        return new Complex(Math.cos(a), Math.sin(a))
    }

    // String representation of a complex number:
    def show() {
        print('('); print(re); print(','); print(im); println(')')
    }
}


/**************************************************************************************************
 * Auxiliary Function: complexShow
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 **************************************************************************************************/
def complexShow(x) {
    N = x.length
    for (k in 0..<N) {
        x[k].show()
    }
}


/**************************************************************************************************
 * Auxiliary Function: timeIt
 *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
 *
 * Parameters:
 *   f
 *     A closure to the Fourier Transform function to be called, receiving one argument and
 *     returning a vector of FT coefficients;
 *   size
 *     Power of two of the size of the vector on which the transform will be applied;
 *   repeat
 *     Number of times the function will be called. Defaults to REPEAT.
 *
 * Returns:
 *   The average execution time for that function with a vector of the given size.
 **************************************************************************************************/
def timeIt(f, size, repeat=REPEAT) {
    x = new Complex[size];                                     // Generate a vector;
    for (j in 0..<size)
        x[j] = new Complex(j as Float, 0.0f);

    t0 = System.currentTimeMillis()                            // Start a timer;
    for (j in 1..repeat)                                       // Repeated calls;
        directFT(x)
    return (System.currentTimeMillis() - t0) / (1000.0*repeat) // Compute average;
}


/**************************************************************************************************
 * Function: directFT
 *   Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
 *   complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the DFT will be computed. Given the nature of the implementation, there
 *     is no restriction on the size of the vector, although it will almost always be called with a
 *     power of two size to give a fair comparison;
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
def directFT(x) {
    N = x.length                                       // Length of the vector;
    X = new Complex[N]                                 // Accumulate the results;

    // Initialize twiddle factors:
    W = Complex.exp(-2*Math.PI/N)
    Wk = new Complex(1.0f, 0.0f)

    for (k in 0..<N) {                                 // Compute the kth coefficient;
        X[k] = new Complex(0.0f, 0.0f)                 // Accumulate the results;
        Wkn = new Complex(1.0f, 0.0f)
        for (n in 0..<N) {                             //   Operate the summation;
            X[k] = X[k] + Wkn * x[n]                   //     Compute every term;
            Wkn = Wkn * Wk                             // Update twiddle factor;
        }
        Wk = Wk * W
    }
    return X
}


/**************************************************************************************************
 * Function: recursiveFFT
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 *   complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector
 *     of a power of two length, or it will fail. No checks on this are made.
 *
 *  Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
def recursiveFFT(x) {
    def N = x.length
    if (N == 1)                                        // A length-1 vector is its own FT;
        return x
    else {
        def N2 = N/2 as Integer
        def xe = new Complex[N2]                       // Allocate memory for computation;
        def xo = new Complex[N2]
        def X = new Complex[N]
        for (k in 0..<N)
            X[k] = new Complex(0.0f, 0.0f)

        for (k in 0..<N2) {                            // Split even and odd samples;
            xe[k] = x[2*k]
            xo[k] = x[2*k+1]
        }
        def Xe = recursiveFFT(xe)                      // Transform of even samples;
        def Xo = recursiveFFT(xo)                      // Transform of odd samples;

        // Twiddle factors:
        W = Complex.exp(-2*Math.PI/N)
        Wk = new Complex(1.0f, 0.0f)
        for (k in 0..<N2) {
            w = Wk * Xo[k]                             // Recombine results;
            X[k] = Xe[k] + w
            X[k+N2] = Xe[k] - w
            Wk = Wk * W                                // Update twiddle factors;
        }
        return X
    }
}


/**************************************************************************************************
 * Method: bitReverse
 *   Bit-reversed version of an integer number.
 *
 * Parameters:
 *   k
 *     The number to be bit-reversed;
 *   r
 *     The number of bits to take into consideration when reversing.
 *
 * Returns:
 *   The number k, bit-reversed according to integers with r bits.
 **************************************************************************************************/
def bitReverse(k, r) {
    def l = 0                                          // Accumulates the results;
    for (i in 0..<r) {                                 // Loop on every bit;
        l = (l << 1) + (k & 1)                         // Tests less signficant bit and add;
        k = k >> 1                                     // Tests next bit;
    }
    return l
}


/**************************************************************************************************
 * Function: iterativeFFT
 *   Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
 *   O(N log_2(N)) complexity, and since there are less function calls, it will probably be
 *   marginally faster than the recursive versions.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector
 *     of a power of two length, or it will fail. No checks on this are made.
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
def iterativeFFT(x) {
    def N = x.length;
    def X = new Complex[N]

    r = Math.round(Math.log2(N)) as Integer                    // Number of bits;
    for (k in 0..<N) {
        l = bitReverse(k, r)                                   // Reorder the vector according to
        X[l] = x[k]                                            //   the bit-reversed order;
    }

    step = 1                                                   // Computation of twiddle factors;
    for (k in 0..<r) {
        l = 0
        while (l < N) {
            W = Complex.exp(-Math.PI/step)                     // Twiddle factors:
            Wkn = new Complex(1.0f, 0.0f)
            for (n in 0..<step) {
                p = l + n
                q = p + step
                X[q] = X[p] - Wkn * X[q]                       // Recombine results;
                X[p] = new Complex(2.0f, 0.0f)*X[p] - X[q]
                Wkn = Wkn * W                                  // Update twiddle factors;
            }
            l = l + 2*step
        }
        step = step << 1
    }
    return X
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
def main() {

    // Start by printing the table with time comparisons:
    println("+---------"*6 + "+")
    println("|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  |")
    println("+---------"*6 + "+")

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for (r in 5..10) {

        // Compute the average execution time:
        n = Math.pow(2, r) as Integer
        dtime = timeIt({ x -> directFT(x) }, n, REPEAT)
        rtime = timeIt({ x -> recursiveFFT(x) }, n, REPEAT)
        itime = timeIt({ x -> iterativeFFT(x) }, n, REPEAT)

        // Print the results:
        println(String.format('| %1$7s | %2$7s | %3$7s | %4$7.7s | %5$7.7s | %6$7.7s |', n, n**2, n*r, dtime, rtime, itime))

    }
    println("+---------"*6 + "+")
}

main()
