/**************************************************************************************************
 * Fast Fourier Transform -- Scala Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 Scala is a compiled language, so the first step is to compile the program. This can be done by
 issuing the command:

 $ scalac fft.scala

 This will generate an executable file in the same folder, that can be run with the command:

 $ scala fft
 **************************************************************************************************/


/**************************************************************************************************
 * Mini-library to deal with complex numbers. While it may be useful to use modules in bigger
 * projects, this one is very small and thorougly tested, and it is simples to just add the
 * code to the main program.
 **************************************************************************************************/
class Complex(val r: Double = 0.0, val i: Double = 0.0) {

    // Add the argument to this, giving the result as a new complex number:
    def +(c: Complex) = new Complex(this.r+c.r, this.i+c.i)

    // Subtract the argument from this, giving the result as a new complex number:
    def -(c: Complex) = new Complex(this.r-c.r, this.i-c.i)

    // Multiply the argument with this, giving the result as a new complex number:
    def *(c: Complex) = new Complex(this.r*c.r - this.i*c.i, this.r*c.i + this.i*c.r)

    // Divide this by the argument, giving the result as a new complex number:
    def /(a: Double) = new Complex(this.r/a, this.i/a)

    // String representation for printing:
    override def toString(): String = f"(${this.r}%7.4f, ${this.i}%7.4f)"
}


/**************************************************************************************************
 * A companion Complex object will allow the instantiation without using `new` keyword, and also
 * compute complex exponentials without instantiation of the Complex class.
 **************************************************************************************************/
object Complex {

    // Initializes the object:
    def apply(r: Double=0.0, i: Double=0.0): Complex = {
        var z = new Complex(r, i)
        z
    }

    // Complex exponential of an angle:
    def exp(a: Double): Complex = {
        return new Complex(Math.cos(a), Math.sin(a))
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
object fft {

    val REPEAT = 500                           // Number of executions to compute average time;

    /**********************************************************************************************
     * Auxiliary Method: complexShow
     *   Pretty printing of an array of complex numbers, used to inspect results.
     *
     * Parameters:
     *   x
     *     A vector of complex numbers, according to the definition above;
     **********************************************************************************************/
     def complexShow(x: Array[Complex])
     {
         for (i <- 0 until x.length)
             println(x(i))
     }


    /**********************************************************************************************
     * Auxiliary Method: timeIt
     *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
     *
     * Parameters:
     *  f
     *    Function to be called, with the given prototype. The first complex vector is the input
     *    vector, the second complex vector is the result of the computation;
     *  size
     *    Number of elements in the vector on which the transform will be applied;
     *  repeats
     *    Number of times the function will be called.
     *
     * Returns:
     *   The average execution time for that function with a vector of the given size.
     **********************************************************************************************/
    def timeIt(f: Array[Complex] => Array[Complex], size: Int, repeats: Int = REPEAT): Double =
    {
        val x = new Array[Complex](size)                       // Initialize the vector;
        for (j <- 0 until size)
            x(j) = Complex(j, 0)

        // Check from here on:
        val t0 = System.nanoTime()                             // Start a timer;
        for (j <- 0 until repeats) {                           // Repeated calls;
            var X = f(x)
            X
        }
        (System.nanoTime() - t0) / (repeats*1.0e9)
    }


    /**********************************************************************************************
     * Method: directFT
     *   Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
     *   complexity.
     *
     * Parameters:
     *   x
     *     The vector of which the DFT will be computed. Given the nature of the implementation,
     *     there is no restriction on the size of the vector, although it will almost always be
     *     called with a power of two size to give a fair comparison;
     *
     * Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    def directFT(x: Array[Complex]): Array[Complex] =
    {
        val N = x.length                               // Length of the vector;
        val X = new Array[Complex](N)                  // Accumulate the results;

        val W = Complex.exp(-2*Math.PI/N)              // Initializes twiddle factors;
        var Wk = Complex(1, 0)

        for (k <- 0 until N) {                         // Compute the kth coefficient;
            X(k) = Complex()                           // Accumulate the results;
            var Wkn = Complex(1, 0)
            for (n <- 0 until N) {                     //   Operate the summation;
                X(k) = X(k) + Wkn * x(n)               //     Compute every term;
                Wkn = Wkn * Wk                         // Update twiddle factor;
            }
            Wk = Wk * W
        }
        X
    }


    /**********************************************************************************************
     * This is a small and very readable DFT implementation using for-expressions. This capability
     * of Scala resembles list comprehensions found in Python and Haskell, and make for very
     * readable, clean and easy interpretable functions. A detailed explanation of every term in the
     * expression follows:
     *
     * ( for (k <- 0 until x.length) yield             // The kth component of the transform;
     *      (for (n <- 0 until x.length) yield         // For every k, the nth term of summation;
     *            x(n) * exp(-2*Math.PI*k*n/x.length)  // The general term of the summation;
     *       ).reduce(_+_) )                           // Sum;
     *       .toArray                                  // Converts to an array;
     * Here is the analysis equation (in LaTeX form) for comparison:
     *
     *  X[k] = \sum_{0}^{N-1} x[n] e^{-j 2 \pi k n / N} \; k = 0 \ldots N-1
     *
     * This implementation is not, of course, very efficient, since it doesn't take advantage of the
     * regularity of the twiddle factors, and computes complex exponentials for every term. Also,
     * it is not as clear as the Python implementation, but it can probably be cleaned up.
     **********************************************************************************************/
    def forFT(x: Array[Complex]): Array[Complex] =
        ( for (k <- 0 until x.length) yield
              (for (n <- 0 until x.length) yield
                   x(n) * Complex.exp(-2*Math.PI*k*n/x.length)).reduce(_+_) ).toArray


    /**********************************************************************************************
     * Method: recursiveFFT
     *   Fast Fourier Transform using a recursive decimation in time algorithm. This has
     *   O(N log_2(N)) complexity.
     *
     * Parameters:
     *   x
     *     The vector of which the FFT will be computed. This should always be called with a vector
     *     of a power of two length, or it will fail. No checks on this are made.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    def recursiveFFT(x: Array[Complex]): Array[Complex] =
    {
        val N = x.length

        if (N==1)                                              // A length-1 vector is its own FT;
            return x
        else {
            val N2 = N >> 1

            // Split even and odd samples;
            val xe = (for (k <- 0 until N2) yield x(2*k)).toArray
            val xo = (for (k <- 0 until N2) yield x(2*k+1)).toArray
            val Xe = recursiveFFT(xe)                          // Transform of even samples;
            val Xo = recursiveFFT(xo)                          // Transform of odd samples;

            val X = new Array[Complex](N)                      // Allocate memory for computation;
            val W = Complex.exp(-2*Math.PI/N)                  // Twiddle factors:
            var Wk = Complex(1, 0)
            for (k <- 0 until N2) {
                X(k) = Xe(k) + Wk*Xo(k)                        // Recombine results;
                X(k+N2) = Xe(k) - Wk*Xo(k)
                Wk = Wk * W                                    // Update twiddle factors;
            }
            X
        }
    }


    /**********************************************************************************************
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
     **********************************************************************************************/
    def bitReverse(k: Int, r: Int): Int =
    {
        var kr = k
        var l = 0                                      // Accumulate the results;
        for (i <- 0 until r) {                         // Loop on every bit;
            l = (l << 1) + (kr & 1)                    // Test less signficant bit and add;
            kr = (kr >> 1)                             // Test next bit;
        }
        return l
    }


    /**********************************************************************************************
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
     **********************************************************************************************/
    def iterativeFFT(x: Array[Complex]): Array[Complex] =
    {
        val N = x.length
        val X = new Array[Complex](N)

        val r = Math.round(Math.log(N)/Math.log(2)).toInt      // Number of precision bits;
        for (k <- 0 to N-1) {
            var l = bitReverse(k, r)                           // Reorder the vector according to
            X(l) = x(k)                                        //   the bit-reversed order;
        }

        var step = 1                                           // Computation of twiddle factors;
        for (k <- 0 until r) {
            for (l <- 0 until N by 2*step) {
                val W = Complex.exp(-Math.PI/step)             // Twiddle factors;
                var Wkn = Complex(1, 0)
                for (n <- 0 until step) {
                    var p = l + n
                    var q = p + step
                    X(q) = X(p) - Wkn*X(q)                     // Recombine results;
                    X(p) = X(p) + X(p) - X(q)
                    Wkn = Wkn * W                              // Update twiddle factors;
                }
            }
            step <<= 1
        }
        X
    }


    /**********************************************************************************************
     * Main Function:
     **********************************************************************************************/
    def main(args: Array[String]) {

        // Start by printing the table with time comparisons:
        println("+---------+---------+---------+---------+---------+---------+---------+")
        println("|    N    |   N^2   | N logN  | Direct  | ForExpr | Recurs. | Intera. |")
        println("+---------+---------+---------+---------+---------+---------+---------+")

        // Try it with vectors with size ranging from 32 to 1024 samples:
        for (r <- 5 to 10) {

            // Compute the average execution time:
            var n = Math.pow(2, r).toInt
            var dtime = timeIt(directFT, n, REPEAT)
            var ftime = timeIt(forFT, n, REPEAT)
            var rtime = timeIt(recursiveFFT, n, REPEAT)
            var itime = timeIt(iterativeFFT, n, REPEAT)

            // Print the results:
            println(f"| ${n}%7d | ${n*n}%7d | ${n*r}%7d | ${dtime}%7.4f | ${ftime}%7.4f | ${rtime}%7.4f | ${itime}%7.4f |")

        }
        println("+---------+---------+---------+---------+---------+---------+---------+")

    }
}
