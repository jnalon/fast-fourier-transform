/**************************************************************************************************
 * Fast Fourier Transform -- Go Version
 * This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
 *
 * José Alexandre Nalon
 **************************************************************************************************
 * To run this file, just type:
 *
 * $ go run anyfft.go
 *
 * If you want to compile to have an executable file, then build it by issuing the command:
 *
 * $ go build anyfft.go
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
 *   Complex exponential of an angle. Convenience function.
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
func DirectFT(x []complex128) []complex128 {
    N := len(x)
    X := make([]complex128, N)                 // Accumulate the results;
    W := CExp(-2*math.Pi/float64(N))           // Initialize twiddle factors:
    Wk := complex(1, 0)
    for k:=0; k<N; k++ {
        Wkn := complex(1, 0)                   // Initialize twiddle factors;
        for n:=0; n<N; n++ {
            X[k] = X[k] + x[n]*Wkn
            Wkn = Wkn * Wk                     // Update twiddle factor;
        }
        Wk = Wk * W
    }
    return X
}


/**************************************************************************************************
 * Function: Factor
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
func Factor(n int) int {
    rn := int(math.Ceil(math.Sqrt(float64(n))))    // Search up to the square root of the number;
    for i:=2; i<=rn; i++ {
        if n%i == 0 {
            return i                               // If remainder is zero, a factor is found;
        }
    }
    return n
}


/**************************************************************************************************
 * Function: RecursiveFFT
 *   Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller
 *   complexity than the direct FT, though the exact value is difficult to compute.
 *
 * Parameters:
 *   x
 *     The vector of which the FFT will be computed. Its length must be a composite number, or else
 *     the computation will be defered to the direct FT, and there will be no efficiency gain.
 *
 * Returns:
 *   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************/
func RecursiveFFT(x []complex128) []complex128 {
    N := len(x)
    X := make([]complex128, N)
    N1 := Factor(N)                            // Smallest prime factor of length;
    if N1 == N {                               // If the length is prime itself,
        return DirectFT(x)                     //   the transform is given by the direct form;
    } else {
        N2 := N / N1                           // Decompose in two factors, N1 being prime;
        xj := make([]complex128, N2)           // Allocate memory for subsequences
        W := CExp(-2*math.Pi/float64(N))       // Twiddle factor;
        Wj := complex(1, 0)
        for j:=0; j<N1; j++ {                  // Compute every subsequence of size N2;
            for n:=0; n<N2; n++ {
                xj[n] = x[n*N1+j]              // Create the subsequence;
            }
            Xj := RecursiveFFT(xj)             // Compute the DFT of the subsequence;
            Wkj := complex(1, 0)
            for k:=0; k<N; k++ {
                X[k] = X[k] + Xj[k%N2] * Wkj   // Recombine results;
                Wkj = Wkj * Wj                 // Update twiddle factors;
            }
            Wj = Wj * W
        }
        return X
    }
}


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
func main() {

    SIZES := [8]int{ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 };

    // Start by printing the table with time comparisons:
    fmt.Println("+---------+---------+---------+---------+")
    fmt.Println("|    N    |   N^2   | Direct  | Recurs. |")
    fmt.Println("+---------+---------+---------+---------+")

    // Try it with vectors with the given sizes:
    for i:=0; i<8; i++ {

        // Compute the average execution time:
        n := SIZES[i]
        dtime := TimeIt(DirectFT, n, REPEAT)
        rtime := TimeIt(RecursiveFFT, n, REPEAT)

        // Print the results:
        fmt.Printf("| %7d | %7d | %7.4f | %7.4f |\n",
                n, n*n, dtime, rtime)
    }

    fmt.Println("+---------+---------+---------+---------+\n")

}
