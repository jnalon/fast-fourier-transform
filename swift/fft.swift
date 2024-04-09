/**************************************************************************************************
 * Fast Fourier Transform -- Swift Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * This program can be executed with the command:
 *
 * $ swift fft.swift
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
 * Function: recursiveFFT
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 *   complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector of 
 *     a power of two length, or it will fail. No checks on this are made.
 *
 *  Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
func recursiveFFT(_ x: [Complex]) -> [Complex] {
    let N = x.count

    if N == 1 {                                                // A length-1 vector is its own FT.
        return x
    } else {
        let N2 = N >> 1
        var xe: [Complex] = []                                 // Allocate memory for computtion;
        var xo: [Complex] = []
        var X = Array(repeating: Complex(real: 0.0, imag: 0.0), count: N)
        for k in 0..<N {                                       // Split even and odd samples;
            if k%2 == 0 {
                xe.append(x[k])
            } else {
                xo.append(x[k])
            }
        }
        let Xe = recursiveFFT(xe)                              // Transform of even samples;
        let Xo = recursiveFFT(xo)                              // Transform of odd samples;

        let W = cexp(-2.0 * Double.pi / Double(N))
        var Wk = Complex(real: 1.0, imag: 0.0)
        for k in 0..<N2 {
            let WXok = Wk * Xo[k]                              // Recombine results;
            X[k] = Xe[k] + WXok
            X[k+N2] = Xe[k] - WXok
            Wk = Wk * W                                        // Update twiddle factors;
        }
        return X
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
func bit_reverse(_ k: Int, _ r: Int) -> Int {
    var k0 = k
    var l: Int = 0                                             // Accumulate the results;
    for _ in 0..<r {                                           // Loop on every bit;
        l = (l << 1) + (k0 & 1)                                // Test significant bit and add;
        k0 = k0 >> 1                                           // Test next bit;
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
 *     The vector of which the FFT will be computed. This should always be called with a vector of
 *     a power of two length, or it will fail. No checks on this are made.
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
func iterativeFFT(_ x: [Complex]) -> [Complex] {
    let N = x.count
    var X = Array(repeating: Complex(real: 0.0, imag:0.0), count: N)

    let r = Int(log(Double(N))/log(2.0))                       // Number of bits;
    for k in 0..<N {
        let l = bit_reverse(k, r)                              // Reorder the vector according to
        X[l] = x[k]                                            //   the bit-reversed order;
    }

    var step = 1                                               // Computation of twiddle factors;
    for _ in 0..<r {
        for l in stride(from: 0, to: N, by: 2*step) {
            let W = cexp(-Double.pi / Double(step))
            var Wkn = Complex(real: 1.0, imag: 0.0)
            for n in 0..<step {
                let p = l + n
                let q = p + step
                X[q] = X[p] - Wkn * X[q]                       // Recombine results;
                X[p] = 2.0 * X[p] - X[q]
                Wkn = Wkn * W                                  // Update twiddle factors;
            }
        }
        step = step << 1
    }
    return X
}


/**************************************************************************************************
 * Main Program.
 **************************************************************************************************/
let REPEAT = 500                               // Number of executions to compute average time;

// Start by printing the table with time comparisons:
print("+---------+---------+---------+---------+---------+---------+")
print("|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |")
print("+---------+---------+---------+---------+---------+---------+")

// Try it with vectors with size ranging from 32 to 1024 samples:
for r in 5...10  {

    // Compute the average execution time:
    let n = Int(pow(2.0, Double(r)))
    let dtime = timeIt(directFT, n, REPEAT)
    let rtime = timeIt(recursiveFFT, n, REPEAT)
    let itime = timeIt(iterativeFFT, n, REPEAT)

    // Print the results:
    print(String(format: "| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |",
                         n, n*n, n*r, dtime, rtime, itime))

}
print("+---------+---------+---------+---------+---------+---------+")
