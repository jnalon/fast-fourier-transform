/**************************************************************************************************
 * Fast Fourier Transform -- Go Version
 * This version implements Cooley-Tukey algorithm for powers of 2 only.
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * To run this file, just type:
 *
 * $ go run fft.go
 *
 * If you want to compile to have an executable file, then build it by issuing the command:
 *
 * $ go build fft.go
 **************************************************************************************************/

package main


/**************************************************************************************************
 Include necessary libraries:
 **************************************************************************************************/
import (
    "fmt"                              // String and output formatting;
    "math"                             // Math functions;
    "time"                             // Time measuring;
)


/**************************************************************************************************
 Definitions:
 **************************************************************************************************/
const REPEAT = 500                     // Number of executions to compute average time;


/**************************************************************************************************
 * Auxiliary function: CExp
 *   Computes the complex exponential of an angle. Convenience function.
 *
 * Parameters:
 *   a
 *     Angle
 *
 * Returns
 *   A complex number with the complex exponential of the angle.
 **************************************************************************************************/
func CExp(a float64) complex128 {
    return complex(math.Cos(a), math.Sin(a))
}


/**************************************************************************************************
 * Auxiliary function: ComplexShow
 *   Pretty printing of an array of complex numbers, used to inspect results.
 *
 * Parameters:
 *   x
 *     A vector of complex numbers, according to the definition above;
 **************************************************************************************************/
func ComplexShow(x []complex128) {
    for i:=0; i<len(x); i++ {
        fmt.Printf("( %7.4f, %7.4f )\n", real(x[i]), imag(x[i]))
    }
}


/**************************************************************************************************
 * Auxiliary function: TimeIt
 *   This function calls a Fast Fourier Transform function repeatedly a certain number of times,
 *   measure execution time and average it.
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
func TimeIt(f func([]complex128) []complex128, size int, repeat int) float64 {
    x := make([]complex128, size)              // Initialize the vector;
    for j:=0; j<size; j++ {
        x[j] = complex(float64(j), 0)
    }
    t0 := time.Now()                           // Starting time;
    for j:=0; j<repeat; j++ {
        f(x)
    }
    t1 := time.Since(t0)
    return time.Duration.Seconds(t1) / float64(repeat)
}


/**************************************************************************************************
 * Function: DirectFt
 *   Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that has
 *   O(N^2) complexity.
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
func DirectFT(x []complex128) []complex128 {
    N := len(x)
    X := make([]complex128, N)                 // Accumulates the results;
    W := CExp(-2*math.Pi/float64(N))           // Initializes twiddle factors:
    Wk := complex(1, 0)
    for k:=0; k<N; k++ {
        Wkn := complex(1, 0)                   // Initializes twiddle factors;
        for n:=0; n<N; n++ {
            X[k] = X[k] + x[n]*Wkn
            Wkn = Wkn * Wk                     // Update twiddle factor;
        }
        Wk = Wk * W
    }
    return X
}


/**************************************************************************************************
 * Function: RecursiveFFT
 *   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This has
 *   O(N log_2(N)) complexity.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector of
 *     a power of two length, or it will fail. No checks on this are made.
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
func RecursiveFFT(x []complex128) []complex128 {
    N := len(x)
    if N == 1 {                                // A length-1 vector is its own FT;
        return x
    } else {
        N2 := N / 2

        xe := make([]complex128, N2)           // Allocates memory for computation;
        xo := make([]complex128, N2)

        for i:=0; i<N/2; i++ {                 // Split even and odd samples;
            xe[i] = x[2*i]
            xo[i] = x[2*i+1]
        }
        Xe := RecursiveFFT(xe)                 // Transform of even samples;
        Xo := RecursiveFFT(xo)                 // Transform of odd samples;

        X := make([]complex128, N)
        W := CExp(-2*math.Pi/float64(N))       // Twiddle factors;
        Wk := complex(1, 0)
        for k:=0; k<N2; k++ {
            w := Wk * Xo[k]                    // Recombine results;
            X[k] = Xe[k] + w
            X[k+N2] = Xe[k] - w
            Wk = Wk * W                        // Update twiddle factors;
        }

        return X
    }
}



/**************************************************************************************************
 * Function: BitReverse
 *   Computes the bit-reversed function of an integer number.
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
func BitReverse(k int, r int) int {
    l := 0                                     // Accumulates the results;
    for i:=0; i<r; i++ {                       // Loop on every bit;
        l = (l << 1) + (k & 1)                 // Tests less signficant bit and add;
        k = k >> 1                             // Tests next bit;
    }
    return l
}


/**************************************************************************************************
 * Function: IterativeFFT
 *   Computes the Fast Fourier Ttransform using an iterative in-place decimation in time algorithm.
 *   This has O(N log_2(N)) complexity, and since there are less function calls, it will probably
 *   be marginally faster than the recursive versions.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. This should always be called with a vector of
 *     a power of two length, or it will fail. No checks on this are made.
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
func IterativeFFT(x []complex128) []complex128 {
    N := len(x)
    r := int(math.Floor(math.Log(float64(N))/math.Log(2)))     // Number of bits;

    X := make([]complex128, N)                 // Reorder the vector according to the
    for k:=0; k<N; k++ {                       //   bit-reversed order;
        l := BitReverse(k, r)
        X[l] = x[k]
    }

    step := 1                                  // Auxiliary for computation of twiddle factors;
    for k:=0; k<r; k++ {
        for l:=0; l<N; l=l+2*step {
            W := CExp(-math.Pi/float64(step))  // Twiddle factors;
            Wkn := complex(1, 0)
            for n:=0; n<step; n++ {
                p := l + n
                q := p + step
                X[q] = X[p] - Wkn * X[q]       // Recombine results;
                X[p] = 2*X[p] - X[q]
                Wkn = Wkn * W                  // Update twiddle factors;
             }
        }
        step = step << 1
    }
    return X
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
func main() {

    // Starts by printing the table with time comparisons:
    fmt.Println("+---------+---------+---------+---------+---------+---------+")
    fmt.Println("|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |")
    fmt.Println("+---------+---------+---------+---------+---------+---------+")

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for r:=5; r<11; r++ {

        // Computes the average execution time:
        n := int(math.Exp2(float64(r)))
        dtime := TimeIt(DirectFT, n, REPEAT)
        rtime := TimeIt(RecursiveFFT, n, REPEAT)
        itime := TimeIt(IterativeFFT, n, REPEAT)

        // Print the results:
        fmt.Printf("| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n",
                n, n*n, r*n, dtime, rtime, itime)
    }

    fmt.Println("+---------+---------+---------+---------+---------+---------+\n")

}
