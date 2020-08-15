/**************************************************************************************************
 * Fast Fourier Transform -- Groovy Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
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
 * Function: factor
 *   Smallest prime factor of a given number. If the argument is prime itself, then it is the
 *   return value.
 *
 * Parameters:
 *   n
 *     Number to be inspected.
 *
 * Returns:
 *   The smallest prime factor, or the number itself if it is already a prime.
 **************************************************************************************************/
def factor(n) {
    rn = n/2;                                  // Search up to half the number;
    for (i in 2..rn)
        if (n%i == 0) return i;                // If remainder is zero, a factor is found;
    return n;
}


/**************************************************************************************************
 * Function: recursiveFFT
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
 *   complexity than the direct FT, though the exact value is difficult to compute.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. Its length must be a composite number, or else
 *     the computation will be defered to the direct FT, and there will be no efficiency gain.
 *
 *  Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
def recursiveFFT(x) {
    def N = x.length
    def N1 = factor(N)                                 // Smallest prime factor of length;
    if (N1 == N)                                       // If the length is prime itself,
        return directFT(x)                             //   transform is given by the direct form;
    else {
        def N2 = (N/N1) as Integer                     // Decompose in two factors, N1 being prime;
        def xj = new Complex[N2]                       // Allocate memory for subsequences;
        def X = new Complex[N]
        for (n in 0..<N)
            X[n] = new Complex(0.0f, 0.0f)

        // Twiddle factors:
        def W = Complex.exp(-2*Math.PI/N)
        def Wj = new Complex(1.0f, 0.0f)
        for (j in 0..<N1) {                            // Compute every subsequence of size N2;
            for (n in 0..<N2)
                xj[n] = x[n*N1+j]                      // Create the subsequence;
            def Xj = recursiveFFT(xj)                  // Compute the DFT of the subsequence;
            Wkj = new Complex(1.0f, 0.0f)
            for (k in 0..<N) {
                X[k] = X[k] + Wkj*Xj[k%N2]             // Recombine results;
                Wkj = Wkj * Wj                         // Update twiddle factors;
            }
            Wj = Wj * W
        }
        return X
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
def main() {
    def SIZES = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ] as Integer[]

    // Start by printing the table with time comparisons:
    println("+---------"*4 + "+")
    println("|    N    |   N^2   | Direct  | Recurs. |")
    println("+---------"*4 + "+")

    // Try it with vectors with the given sizes:
    for (i in 0..<SIZES.length) {

        // Compute the average execution time:
        n = SIZES[i]
        dtime = timeIt({ x -> directFT(x) }, n, REPEAT)
        rtime = timeIt({ x -> recursiveFFT(x) }, n, REPEAT)

        // Print the results:
        println(String.format('| %1$7s | %2$7s | %3$7.7s | %4$7.7s |', n, n**2, dtime, rtime))

    }
    println("+---------"*6 + "+")
}

main()
