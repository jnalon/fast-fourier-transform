/**************************************************************************************************
 * Fast Fourier Transform -- Dart Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program can be run by issuing the command:
 *
 * $ dart fft.dart
 *
 **************************************************************************************************/

/**************************************************************************************************
 * Include necessary libraries:
 **************************************************************************************************/
import 'dart:math';                            // Math functions;


/**************************************************************************************************
 * Mini-library to deal with complex numbers.
 **************************************************************************************************/
class Complex {

    double r;
    double i;

    // Constructor:
    Complex(this.r, this.i);

    // Add the argument to this, giving the result as a new complex number:
    Complex operator +(Complex c) => new Complex(r + c.r, i + c.i);

    // Subtract the argument from this, giving the result as a new complex number:
    Complex operator -(Complex c) => new Complex(r - c.r, i - c.i);

    // Multiply the argument with this, giving the result as a new complex number:
    Complex operator *(Complex c) => new Complex(r*c.r - i*c.i, r*c.i + i*c.r);

    // Divide this by the argument, giving the result as a new complex number:
    Complex operator /(double a) => new Complex(r/a, i/a);

    // Convert to string:
    String toString() => "(${r}, ${i})";

}

// Complex exponential of an angle:
Complex Cexp(double a) => new Complex(cos(a), sin(a));


/**************************************************************************************************
 * Auxiliary functions to format numbers for printing on screen.
 * Both receive the numbers to be formatted and output a string ready to be printed.
 **************************************************************************************************/
String fi(int d)
{
    return d.toString().padLeft(7, ' ');
}

String fd(double f)
{
    f = (10000 * f).round() / 10000;
    var l = f.toString().split('.');
    return l[0].padLeft(2, ' ') + '.' + l[1].padRight(4, '0').substring(0, 4);
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
double timeIt(Function f, int size, int repeat)
{
    var x = List<Complex>.generate(size, (int i) => new Complex(i.toDouble(), 0.0));
    final stopwatch = Stopwatch()..start();
    for (var j=0; j<repeat; j++)
        f(x);
    stopwatch.stop();
    return stopwatch.elapsedMilliseconds / (1000.0 * repeat);
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
List<Complex> directFT(List<Complex> x)
{
    var N = x.length;
    var X = List<Complex>.filled(N, Complex(0.0, 0.0));        // Accumulate the results;
    var W = Cexp(-2*pi/N);                                     // Initialize twiddle factors;
    var Wk = Complex(1.0, 0.0);
    for (var k=0; k<N; k++) {
        var Wkn = Complex(1.0, 0.0);
        for (var n=0; n<N; n++) {
            X[k] = X[k] + Wkn * x[n];
            Wkn = Wkn * Wk;                                    // Update twiddle factor;
        }
        Wk = Wk * W;
    }
    return X;                                                  // Return value;
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
List<Complex> recursiveFFT(List<Complex> x)
{
    var N = x.length;
    if (N==1) {                                                // A length-1 vector is its own FT;
        return x;
    } else {
        int N2 = (N / 2).toInt();
        var X = List<Complex>.filled(N, Complex(0.0, 0.0));    // Allocate memory for computation;
        var xe = List<Complex>.generate(N2, (int n) => x[2*n]);
        var xo = List<Complex>.generate(N2, (int n) => x[2*n+1]);
        var Xe = recursiveFFT(xe);                             // Transform of even samples;
        var Xo = recursiveFFT(xo);                             // Transform of odd samples;

        var W = Cexp(-2*pi/N);                                 // Twiddle factors;
        var Wk = Complex(1.0, 0.0);
        for (var k=0;  k<N2; k++) {
            var w = Wk * Xo[k];                                // Recombine results;
            X[k] = Xe[k] + w;
            X[k+N2] = Xe[k] - w;
            Wk = Wk * W;                                       // Update twiddle factors;
        }
        return X;                                              // Return value;
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
int bitReverse(int k, int r)
{
    var l = 0;                                         // Accumulate the results;
    for (var i=0; i<r; i++) {                          // Loop on every bit;
        l = (l << 1) + (k % 2);                        // Test less signficant bit and add;
        k = (k >> 1);                                  // Test next bit;
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
List<Complex> iterativeFFT(List<Complex> x)
{
    var N = x.length;
    var X = List<Complex>.filled(N, Complex(0, 0));
    var r = (log(N)/ln2).round();                      // Number of bits;
    for (var k=0; k<N; k++) {
        var l = bitReverse(k, r);                      // Reorder the vector according to
        X[l] = x[k];                                   //   the bit-reversed order;
    }

    var step = 1;                                      // Computation of twiddle factors;
    for (var k=0; k<r; k++) {
        var l = 0;
        while (l < N) {
            var W = Cexp(-pi/step);                    // Twiddle factors;
            var Wkn = Complex(1.0, 0.0);
            for (var n=0; n<step; n++) {
                var p = l + n;
                var q = p + step;
                X[q] = X[p] - Wkn * X[q];              // Recombine results;
                X[p] = X[p]*Complex(2.0, 0.0) - X[q];
                Wkn = Wkn * W;                         // Update twiddle factors;
            }
            l = l + 2*step;
        }
        step = step * 2;
    }

    return X;                                          // Return value;
}



/**************************************************************************************************
 * Main function:
 **************************************************************************************************/
void main() {

    const int repeat = 500;                    // Number of executions to compute average time;

    // Start by printing the table with time comparisons:
    print("+---------+---------+---------+---------+---------+---------+");
    print("|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |");
    print("+---------+---------+---------+---------+---------+---------+");

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for (var r=5; r<=10; r++) {

        // Compute the average execution time:
        var n = pow(2, r).toInt();
        var dtime = timeIt(directFT, n, repeat);
        var rtime = timeIt(recursiveFFT, n, repeat);
        var itime = timeIt(iterativeFFT, n, repeat);

        // Print the results:
        var results = ("| ${fi(n)} | ${fi(n*n)} | ${fi(n*r)} | ${fd(dtime)} | ${fd(rtime)} | ${fd(itime)} |");
        print(results);
    }

    print("+---------+---------+---------+---------+---------+---------+");
}
