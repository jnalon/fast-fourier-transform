/**************************************************************************************************
 * Fast Fourier Transform -- Dart Version
 * This version implements Cooley-Tukey algorithm for composite-length sequences.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program can be run by issuing the command:
 *
 * $ dart anyfft.dart
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
int factor(int n)
{
    var rn = (sqrt(n)).ceil();                 // Search up to the square root of the number;
    for (var i=2; i<=rn; i++)
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
List<Complex> recursiveFFT(List<Complex> x)
{
    var N = x.length;
    var N1 = factor(N);                                        // Smallest prime factor of length;
    if (N == N1) {                                             // If the length is prime itself,
        return directFT(x);                                    //   transform is given by the direct form;
    } else {
        var N2 = (N / N1).toInt();                             // Decompose in two factors, N1 being prime;
        var X = List<Complex>.filled(N, Complex(0.0, 0.0));    // Allocate memory for computation;
        var W = Cexp(-2*pi/N);                                 // Twiddle factor;
        var Wj = Complex(1.0, 0.0);
        for (var j=0; j<N1; j++) {                             // Compute every subsequence of size N2;
            var xj = List<Complex>.generate(N2, (int n) => x[n*N1 + j]);  // Create the subsequence;
            var Xj = recursiveFFT(xj);                         // Compute the DFT of the subsequence;
            var Wkj = Complex(1.0, 0.0);
            for (var k=0; k<N; k++) {
                X[k] = X[k] + Xj[k%N2] * Wkj;                  // Recombine results;
                Wkj = Wkj * Wj;                                // Update twiddle factors;
            }
            Wj = Wj * W;
        }
        return X;
    }
}


void main() {

    var sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ];
    const int repeat = 500;                    // Number of executions to compute average time;

    // Start by printing the table with time comparisons:
    print("+---------+---------+---------+---------+");
    print("|    N    |   N^2   | Direct  | Recurs. |");
    print("+---------+---------+---------+---------+");

    // Try it with vectors with the given sizes:
    for (var n in sizes) {

        // Compute the average execution time:
        var dtime = timeIt(directFT, n, repeat);
        var rtime = timeIt(recursiveFFT, n, repeat);

        // Print the results:
        var results = "| ${fi(n)} | ${fi(n*n)} | ${fd(dtime)} | ${fd(rtime)} |";
        print(results);
    }

    print("+---------+---------+---------+---------+");

}
