/**************************************************************************************************
 * Fast Fourier Transform -- Scala Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 Scala is a compiled language, so the first step is to compile the program. This can be done by
 issuing the command:

 $ scalac anyfft.scala

 This will generate an executable file in the same folder, that can be run with the command:

 $ scala anyfft
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
object anyfft {

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
     * Method: factor
     *   Smallest prime factor of a given number. If the argument is prime itself, then it is the
     *   return value.
     *
     * Parameters:
     *   n
     *     Number to be inspected.
     *
     * Returns:
     *   The smallest prime factor, or the number itself if it is already a prime.
     **********************************************************************************************/
    def factor(n: Int): Int =
    {
        val rn = n/2                                   // Search up to half the number;
        for (i <- 2 to rn)
            if (n%i == 0) return i                     // If remainder is zero, a factor is found;
        n
    }


    /**********************************************************************************************
     * Method: recursiveFFT
     *   Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
     *   complexity than the direct FT, though the exact value is difficult to compute.
     *
     * Parameters:
     *   x
     *   The vector of which the FFT will be computed. Its length must be a composite number, or
     *   else the computation will be defered to the direct FT, and there will be no efficiency
     *   gain.
     *
     *  Returns:
     *   A complex-number vector of the same size, with the coefficients of the DFT.
     **********************************************************************************************/
    def recursiveFFT(x: Array[Complex]): Array[Complex] =
    {
        val N = x.length
        val N1 = factor(N)                             // Smallest prime factor of length;

        if (N1==N)                                     // If the length is prime itself,
            return directFT(x)                         //   transform is given by the direct form;
        else {
            val N2 = N / N1

            val X = new Array[Complex](N)              // Allocate memory for computation;
            for (n <- 0 until N)
                X(n) = new Complex()

            val W = Complex.exp(-2*Math.PI/N)          // Twiddle factors;
            var Wj = new Complex(1, 0)
            for (j <- 0 until N1) {                    // Compute every subsequence of size N2;
                val xj = (for (n <- 0 until N2) yield x(n*N1+j)).toArray
                val Xj = recursiveFFT(xj)              // Compute the DFT of the subsequence;
                var Wkj = new Complex(1, 0)
                for (k <- 0 until N) {
                    X(k) = X(k) + Wkj*Xj(k%N2)         // Recombine results;
                    Wkj = Wkj * Wj                     // Update twiddle factors;
                }
                Wj = Wj * W
            }
            X
        }
    }


    /**********************************************************************************************
     * Main Function:
     **********************************************************************************************/
    def main(args: Array[String]) {

        val SIZES = Array( 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 )

        // Start by printing the table with time comparisons:
        println("+---------+---------+---------+---------+")
        println("|    N    |   N^2   | Direta  | Recurs. |")
        println("+---------+---------+---------+---------+")

        // Try it with vectors with the given sizes:
        for (n <- SIZES) {

            // Compute the average execution time:
            var dtime = timeIt(directFT, n, REPEAT)
            var rtime = timeIt(recursiveFFT, n, REPEAT)

            // Print the results:
            println(f"| ${n}%7d | ${n*n}%7d | ${dtime}%7.4f | ${rtime}%7.4f |")

        }
        println("+---------+---------+---------+---------+")

    }
}
