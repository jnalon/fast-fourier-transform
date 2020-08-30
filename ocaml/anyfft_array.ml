(**************************************************************************************************
 Fast Fourier Transform -- OCaml Version
 This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only) using
 OCaml arrays.

 José Alexandre Nalon
 **************************************************************************************************
 There are different ways to compile and run an OCaml program. To run directly from the command
 line, just type:

 $ ocaml anyfft_array.ml

 And the program will run. To compile the program, use the command:

 $ ocamlc -o anyfft_array anyfft_array.ml

 This will generate an executable file named `anyfft_array` in the same directory, that can be run
 from the command line by typing:

 $ ./anyfft_array
 **************************************************************************************************)

(**************************************************************************************************
 Open needed modules:
 **************************************************************************************************)
open Array                             (* Array of numbers; *)
open Complex                           (* Complex numbers; *)
open Printf                            (* Output results; *)
open Sys                               (* Time measurement; *)


(**************************************************************************************************
 Definitions:
 **************************************************************************************************)
let repeats = 50
let pi = 3.1415926535897931


(**************************************************************************************************
 Auxiliary Functions:
 **************************************************************************************************)
(* Convert integer to complex. *)
let complex_of_int n = { re=(float_of_int n); im=0. }


(* Integer power of a complex number. *)
let cpow w n = Complex.pow w (complex_of_int n)


(* Creates a sequence of integer numbers. *)
let range a b = List.init (b - a + 1) (fun n -> n + a)


(**************************************************************************************************
 Auxiliary function: complex_show
   Pretty printing of an array of complex numbers, used to inspect results.

 Parameters:
   x
     An array of complex numbers;
 **************************************************************************************************)
let complex_show x =
    for i = 0 to ((length x) - 1) do
        printf "( %7.4f, %7.4f )\n" x.(i).re x.(i).im
    done


(**************************************************************************************************
 Measure execution time through repeated calls to a (Fast) Fourier Transform function.

 Parameters:
   f
     Function to be called;
   size
     Power of two of the size of the vector on which the transform will be applied;
   repeats
     Number of times the function will be called.

 Returns:
   The average execution time for that function with a vector of the given size.
 **************************************************************************************************)
let time_it f size repeats =
    let x = Array.init size (fun n -> { re=(float n); im=0. }) in
    let t0 = Sys.time () in
    for i = 1 to size do
        f x;
    done;
    (Sys.time () -. t0) /. (float repeats)


(**************************************************************************************************
 Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.

 Parameters:
   x
     The vector of which the DFT will be computed. Given the nature of the implementation, there is
     no restriction on the size of the list, although it will almost always be called with a power
     of two size to give a fair comparison.

 Returns:
     A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let direct_ft x =
    let nx = length x in
    let tX = Array.init nx (fun n -> Complex.zero) in
    let w = exp { re=0.; im= -2. *. pi /. (float nx) } in
    let wk = ref { re=1.; im=0. } in
    for k = 0 to (nx - 1) do
        let wkn = ref { re=1.; im=0. } in
        for n = 0 to (nx - 1) do
            tX.(k) <- add tX.(k) (mul !wkn x.(n));
            wkn := mul !wkn !wk
        done;
        wk := mul !wk w
    done;
    tX


(**************************************************************************************************
 Smallest prime factor of a given number. If the argument is prime itself, then it is the return
 value.

 Parameters:
   n
     Number to be inspected;

 Returns:
   The smallest prime factor, or the number itself if it is already a prime.
 **************************************************************************************************)
let factor n =
    let rec loop n k =
        if n mod k == 0 then
            k
        else if k == n / 2 then
            n
        else
            loop n (k + 1)
    in loop n 2


(**************************************************************************************************
 Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 complexity.

 Parameters:
   x
     The vector of which the FFT will be computed. Its length must be a composite number, or else
     the computation will be defered to the direct FT, and there will be no efficiency gain.

 Returns:
   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let rec recursive_fft x =
    let nx = Array.length x in
    let n1 = factor nx in
    if n1 == nx then
        direct_ft x
    else
        let n2 = nx / n1 in
        let xj = Array.make n2 Complex.zero in
        let tX = Array.make nx Complex.zero in
        let w = exp { re=0.; im= -2. *. pi /. (float nx) } in
        let wj = ref { re=1.; im=0. } in
        for j = 0 to (n1 - 1) do
            for k = 0 to (n2 - 1) do
                xj.(k) <- x.(k * n1 + j)
            done;
            let tXj = recursive_fft xj in
            let wkj = ref { re=1.; im=0. } in
            for k = 0 to (nx - 1) do
                tX.(k) <- add tX.(k) (mul tXj.(k mod n2) !wkj);
                wkj := mul !wkj !wj
            done;
            wj := mul !wj w
        done;
        tX


(**************************************************************************************************
 Main function:
 **************************************************************************************************)
let main () =

    (* Start by printing the table with time comparisons: *)
    printf "+---------+---------+---------+---------+\n";
    printf "|    N    |   N^2   | Direta  | Recurs. |\n";
    printf "+---------+---------+---------+---------+\n";

    (* Try it with vectors with the given sizes: *)
    let sizes = [ 2*3; 2*2*3; 2*3*3; 2*3*5; 2*2*3*3; 2*2*5*5; 2*3*5*7; 2*2*3*3*5*5 ] in
    let main_loop n =

        (* Compute the average execution time: *)
        let dtime = time_it direct_ft n repeats in
        let rtime = time_it recursive_fft n repeats in

        (* Print the results: *)
        printf "| %7d | %7d | %7.4f | %7.4f |\n" n (n*n) dtime rtime;
        flush stdout

    in List.iter main_loop sizes;

    printf "+---------+---------+---------+---------+\n"

let () = main ()
