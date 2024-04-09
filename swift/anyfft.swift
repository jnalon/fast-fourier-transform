/**************************************************************************************************
 * Fast Fourier Transform -- Swift Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program can be executed with the command:
 *
 * $ swift anyfft.swift
 *
 **************************************************************************************************/

/**************************************************************************************************
 * Include necessary libraries:
 **************************************************************************************************/
import Foundation                                              // Basic functionality;

/**************************************************************************************************
 * Mini-library to deal with complex numbers. While it may be useful to use modules in bigger
 * projects, this one is very small and thorougly tested, and it is simples to just add the
 * code to the main program.
 **************************************************************************************************/
struct Complex {

    var real: Double
    var imag: Double

    // Adds two complex numbers:
    static func + (a: Complex, b: Complex) -> Complex {
        return Complex(real: a.real+b.real, imag: a.imag+b.imag)
    }

    // Subtract one complex number from the other:
    static func - (a: Complex, b: Complex) -> Complex {
        return Complex(real: a.real-b.real, imag: a.imag-b.imag)
    }

    // Multiplies two complex numbers:
    static func * (a: Complex, b: Complex) -> Complex {
        return Complex(real: a.real*b.real - a.imag*b.imag,
                       imag: a.real*b.imag + a.imag*b.real)
    }

    // Multiplies a double-precision scalar with a complex number:
    static func * (a: Double, b: Complex) -> Complex {
        return Complex(real: a * b.real, imag: a * b.imag)
    }

    // Division of a complex number by a double-precision scalar:
    static func / (a: Complex, b: Double) -> Complex {
        return Complex(real: a.real / b, imag: a.imag / b)
    }

}

// Complex exponention of an angle:
func cexp(_ a: Double) -> Complex {
    return Complex(real: cos(a), imag: sin(a))
}


/**************************************************************************************************
 * Auxiliary Function: complexShow
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 **************************************************************************************************/
func complexShow(_ x: [Complex]) {
    for (index, value) in x.enumerated() {
        print("\(index) -> \(value)")
    }
}


/**************************************************************************************************
 * Auxiliary Function: timeIt
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
 **************************************************************************************************/
func timeIt(_ f: ([Complex])->[Complex], _ size: Int, _ repeats: Int) -> Double {
    var x: [Complex] = []                                      // Initialize the vector;
    for i in 0..<size {
        x.append(Complex(real: Double(i), imag: 0))
    }
    let clock = ContinuousClock()
    let ttime = clock.measure {                                // Measure the time;
        for _ in 0..<repeats {                                 // Repeated calls;
            _ = f(x)
        }
    }
    let double_ttime = Double(ttime.components.seconds)
                       + Double(ttime.components.attoseconds) * 1e-18
    return double_ttime / Double(repeats)
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
func directFT(_ x: [Complex]) -> [Complex] {
    let N = x.count                                            // Length of the vector;
    var X = Array(repeating: Complex(real: 0.0, imag: 0.0), count: N)  // Accumulate the results/

    let W = cexp(-2.0 * Double.pi / Double(N))
    var Wk = Complex(real:1.0, imag:0.0)

    for k in 0..<N {                                           // Compute the kth coefficient;
        var Wkn = Complex(real:1.0, imag:0.0)                  // Accumulate the results;
        for n in 0..<N {                                       //   Operate the summation;
            X[k] = X[k] + Wkn * x[n]                           //     Compute every term;
            Wkn = Wkn * Wk                                     // Update twoddle factor;
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
func factor(_ n: Int) -> Int {
    let rn = Int(ceil(sqrt(Double(n))))
    for i in 2...rn {
        if n % i == 0 {
            return i
        }
    }
    return n
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
func recursiveFFT(_ x: [Complex]) -> [Complex] {
    let N = x.count
    let N1 = factor(N)
    if N1 == N {
        return directFT(x)
    } else {
        let N2 = N / N1

        var xj = Array(repeating: Complex(real: 0.0, imag: 0.0), count: N2)
        var X = Array(repeating: Complex(real: 0.0, imag: 0.0), count: N)

        let W = cexp(-2.0 * Double.pi / Double(N))             // Twiddle factors;
        var Wj = Complex(real: 1.0, imag: 0.0)
        for j in 0..<N1 {
            for n in 0..<N2 {
                xj[n] = x[n*N1+j]
            }
            let Xj = recursiveFFT(xj)
            var Wkj = Complex(real: 1.0, imag: 0.0)
            for k in 0..<N {
                X[k] = X[k] + Wkj * Xj[k%N2]
                Wkj = Wkj * Wj
            }
            Wj = Wj * W
        }
        return X
    }
}


/**************************************************************************************************
 * Main Program.
 **************************************************************************************************/
let SIZES = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]
let REPEAT = 500                               // Number of executions to compute average time;

// Start by printing the table with time comparisons:
print("+---------+---------+---------+---------+")
print("|    N    |   N^2   | Direct  | Recurs. |")
print("+---------+---------+---------+---------+")

// Try it with vectors with size ranging from 32 to 1024 samples:
for n in SIZES {

    // Compute the average execution time:
    let dtime = timeIt(directFT, n, REPEAT)
    let rtime = timeIt(recursiveFFT, n, REPEAT)

    // Print the results:
    print(String(format: "| %7d | %7d | %7.4f | %7.4f |",
                         n, n*n, dtime, rtime))

}
print("+---------+---------+---------+---------+")
