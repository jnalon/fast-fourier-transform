/**************************************************************************************************
 * Fast Fourier Transform -- Kotlin Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program can be compiled by issuing the command:
 *
 * $ kotlinc fft.kt
 *
 * It will generate a file named 'FftKt.class' in the same directory. It can be run by issuing the
 * command:
 *
 * $ kotlin FftKt
 *
 **************************************************************************************************/

/**************************************************************************************************
 * Include necessary libraries:
 **************************************************************************************************/
import kotlin.math.*                   // Math functions;


/**************************************************************************************************
 * Mini-library to deal with complex numbers.
 **************************************************************************************************/
class Complex(val r: Double, val i: Double) {

    // Constructor:
    constructor():this(0.0, 0.0) { }

    // Add the argument to this, giving the result as a new complex number:
    operator fun plus(c: Complex):Complex
    {
        return Complex(r + c.r, i + c.i)
    }

    // Subtract the argument from this, giving the result as a new complex number:
    operator fun minus(c: Complex): Complex
    {
        return Complex(r - c.r, i - c.i)
    }

    // Multiply the argument with this, giving the result as a new complex number:
    operator fun times(c: Complex): Complex
    {
        return Complex(r*c.r - i*c.i, r*c.i + i*c.r)
    }

    // Multiply with an scalar, giving the reulst as a new complex number:
    operator fun times(a: Double): Complex
    {
        return Complex(a*r, a*i)
    }

    // Divide this by the argument, giving the result as a new complex number:
    operator fun div(a: Double): Complex
    {
        return Complex(r/a, i/a)
    }

}

// Complex exponential of an angle:
fun Cexp(a: Double): Complex
{
    return Complex(cos(a), sin(a));
}


/**************************************************************************************************
 * Auxiliary Function: complexShow
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 **************************************************************************************************/
fun complexShow(x: Array<Complex>)
{
    for (i in 0..x.size-1)
        println("( " + x[i].r + ", " + x[i].i +")")
}


/**************************************************************************************************
 * Auxiliary Function: timeIt
 *   Measure execution time through repeated calls to a (Fast) Fourier Transform function.
 *
 * Parameters:
 *  f
 *    Function to be called, with the given prototype. The first complex vector is the input
 *    vector, the second complex vector is the result of the computation, and the integer is the
 *    number of elements in the vector;
 *  size
 *    Number of elements in the vector on which the transform will be applied;
 *  repeat
 *    Number of times the function will be called.
 *
 * Returns:
 *   The average execution time for that function with a vector of the given size.
 **************************************************************************************************/
fun timeIt(f: (x: Array<Complex>) -> Array<Complex>, size: Int, repeat: Int): Double
{
    val x = Array<Complex>(size) { i -> Complex(i.toDouble(), 0.0) }
    val start = System.currentTimeMillis()
    for (j in 1..repeat) {
        f(x)
    }
    return (System.currentTimeMillis() - start).toDouble() / (1000*repeat).toDouble()
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
fun directFT(x: Array<Complex>): Array<Complex>
{
    val N = x.size
    val X = Array<Complex>(N) { _ -> Complex() }       // Accumulate the results;

    val W = Cexp(-2*PI/N.toDouble())                   // Initialize twiddle factors;
    var Wk = Complex(1.0, 0.0)

    for (k in 0..N-1) {
        var Wkn = Complex(1.0, 0.0)
        for (n in 0..N-1) {
            X[k] = X[k] + Wkn * x[n]
            Wkn = Wkn * Wk                             // Update twiddle factor;
        }
        Wk = Wk * W
    }
    return X                                           // Return value;
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
fun recursiveFFT(x: Array<Complex>): Array<Complex>
{
    val N = x.size

    if (N==1) {                                        // A length-1 vector is its own FT;
        return x
    } else {
        val N2 = N / 2

        val X = Array<Complex>(N) { _ -> Complex() }   // Allocate memory for computation;
        val xe = Array<Complex>(N/2) { n -> x[2*n] }
        val xo = Array<Complex>(N/2) { n -> x[2*n+1] }
        val Xe = recursiveFFT(xe)                      // Transform of even samples;
        val Xo = recursiveFFT(xo)                      // Transform of odd samples;

        val W = Cexp(-2*PI/N.toDouble())               // Twiddle factors;
        var Wk = Complex(1.0, 0.0)
        for (k in 0..N2-1) {
            val w = Wk * Xo[k]                         // Recombine results;
            X[k] = Xe[k] + w
            X[k+N2] = Xe[k] - w
            Wk = Wk * W                                // Update twiddle factors;
        }
        return X                                       // Return value;
    }
}


/**************************************************************************************************
 * Function: bitReverse
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
fun bitReverse(k: Int, r: Int): Int
{
    var l: Int = 0;                                    // Accumulate the results;
    var k0: Int = k;
    for (i in 1..r) {                                  // Loop on every bit;
        l = (l shl 1) + (k0 and 1);                    // Test less signficant bit and add;
        k0 = (k0 shr 1);                               // Test next bit;
    }
    return l;
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
fun iterativeFFT(x: Array<Complex>): Array<Complex>
{
    val N = x.size
    val X = Array<Complex>(N) { _ -> Complex() }

    val r = (round(log(N.toDouble(), 2.0))).toInt()    // Number of bits;
    for (k in 0..N-1) {
        var l = bitReverse(k, r)                       // Reorder the vector according to
        X[l] = x[k]                                    //   the bit-reversed order;
    }

    var step: Int = 1                                  // Computation of twiddle factors;
    for (k in 1..r) {
        var l = 0
        while (l < N) {
            var W = Cexp(-PI/step.toDouble())          // Twiddle factors;
            var Wkn = Complex(1.0, 0.0)
            for (n in 0..step-1) {
                var p = l + n
                var q = p + step
                X[q] = X[p] - Wkn * X[q]               // Recombine results;
                X[p] = X[p]*2.0 - X[q]
                Wkn = Wkn * W                          // Update twiddle factors;
            }
            l = l + 2*step
        }
        step = step * 2
    }

    return X                                           // Return value;
}


/**************************************************************************************************
 * Main function:
 **************************************************************************************************/
fun main() {

    val repeat: Int = 500                      // Number of executions to compute average time;

    // Start by printing the table with time comparisons:
    println("+---------+---------+---------+---------+---------+---------+")
    println("|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |")
    println("+---------+---------+---------+---------+---------+---------+")

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for (r in 5..10) {

        // Compute the average execution time:
        var n: Int = 2.0.pow(r).toInt()
        var dtime = timeIt(::directFT, n, repeat)
        var rtime = timeIt(::recursiveFFT, n, repeat);
        var itime = timeIt(::iterativeFFT, n, repeat);

        // Print the results:
        val results = "| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |".format(n, n*n, n*r, dtime, rtime, itime);
        println(results);
    }

    println("+---------+---------+---------+---------+---------+---------+");

}

