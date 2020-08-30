(**************************************************************************************************
 Fast Fourier Transform -- OCaml Version
 This version implements Cooley-Tukey algorithm for powers of 2 only using OCaml arrays.

 José Alexandre Nalon
 **************************************************************************************************
 There are different ways to compile and run an OCaml program. To run directly from the command
 line, just type:

 $ ocaml fft_array.ml

 And the program will run. To compile the program, use the command:

 $ ocamlc -o fft_array fft_array.ml

 This will generate an executable file named `fft_array` in the same directory, that can be run
 from the command line by typing:

 $ ./fft_array
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
 Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 complexity.

 Parameters:
   x
     The vector of which the DFT will be computed. This should have length power of 2;

 Returns:
   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let rec recursive_fft x =
    if Array.length x == 1 then
        x
    else
        let nx = length x in
        let n2 = nx / 2 in
        let xe = Array.make n2 Complex.zero in
        let xo = Array.make n2 Complex.zero in
        let tX = Array.make nx Complex.zero in
        for k = 0 to (n2 - 1) do
            xe.(k) <- x.(2*k);
            xo.(k) <- x.(2*k + 1)
        done;
        let tXe = recursive_fft xe in
        let tXo = recursive_fft xo in
        let w = exp { re=0.; im= -2. *. pi /. (float nx) } in
        let wk = ref { re=1.; im=0. } in
        for k = 0 to (n2 - 1) do
            let wkxo = mul !wk tXo.(k) in
            tX.(k) <- add tXe.(k) wkxo;
            tX.(k+n2) <- sub tXe.(k) wkxo;
            wk := mul !wk w
        done;
        tX


(**************************************************************************************************
 Bit-reversed version of an integer number.

 Parameters:
   k
     The number to be bit-reversed;
   r
     The number of bits to take into consideration when reversing.

 Returns:
   The number k, bit-reversed according to integers with r bits.
 **************************************************************************************************)
let bit_reverse k r =
    let l = ref 0 in
    let k = ref k in
    for i = 1 to r do
        l := 2 * !l + (!k mod 2);
        k := !k / 2
    done;
    !l


(**************************************************************************************************
 Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
 O(N log_2(N)) complexity, and since there are less function calls, it will probably be marginally
 faster than the recursive versions.

 Parameters:
   x
     The vector of which the FFT will be computed. This should always be called with a vector of a
     power of two length, or it will fail. No checks on this are made.

 Returns:
   A complex-number vector of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let iterative_fft x =
    let nx = length x in
    let tX = Array.make nx Complex.zero in
    let r = int_of_float ((Float.log (float nx)) /. (Float.log 2.)) in
    for k = 0 to (nx - 1) do
        let l = bit_reverse k r in
        tX.(l) <- x.(k);
    done;
    let step = ref 1 in
    for k = 0 to (r - 1) do
        let l = ref 0 in
        while !l < nx do
            let w = exp { re=0.; im= -. pi /. (float !step) } in
            let wkn = ref { re=1.; im=0. } in
            for n = 0 to (!step - 1) do
                let p = !l + n in
                let q = p + !step in
                tX.(q) <- sub tX.(p) (mul !wkn tX.(q));
                tX.(p) <- sub (mul { re=2.; im=0. } tX.(p)) tX.(q);
                wkn := mul !wkn w
            done;
            l := !l + 2 * !step
        done;
        step := !step * 2
    done;
    tX


(**************************************************************************************************
 Main function:
 **************************************************************************************************)
let main () =

    (* Start by printing the table with time comparisons: *)
    printf "+---------+---------+---------+---------+---------+---------+\n";
    printf "|    N    |   N^2   | N logN  | Direta  | Recurs. | Itera.  |\n";
    printf "+---------+---------+---------+---------+---------+---------+\n";

    (* Try it with vectors with size ranging from 32 to 1024 samples: *)
    for r = 5 to 10 do

        (* Compute the average execution time: *)
        let n = int_of_float (2. ** (float r)) in
        let dtime = time_it direct_ft n repeats in
        let rtime = time_it recursive_fft n repeats in
        let itime = time_it iterative_fft n repeats in

        (* Print the results: *)
        printf "| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n" n (n*n) (n*r) dtime rtime itime;
        flush stdout

    done;
    printf "+---------+---------+---------+---------+---------+---------+\n"


let () = main ()
