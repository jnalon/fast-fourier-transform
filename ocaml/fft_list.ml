(**************************************************************************************************
 Fast Fourier Transform -- OCaml Version
 This version implements Cooley-Tukey algorithm for powers of 2 only using OCaml lists.

 José Alexandre Nalon
 **************************************************************************************************
 There are different ways to compile and run an OCaml program. To run directly from the command
 line, just type:

 $ ocaml fft_list.ml

 And the program will run. To compile the program, use the command:

 $ ocamlc -o fft_list fft_list.ml

 This will generate an executable file named `fft_list` in the same directory, that can be run from
 the command line by typing:

 $ ./fft_list
 **************************************************************************************************)

(**************************************************************************************************
 Open needed modules:
 **************************************************************************************************)
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
   Pretty printing of a list of complex numbers, used to inspect results.

 Parameters:
   x
     A list of complex numbers;
 **************************************************************************************************)
let rec complex_show x =
    match x with
    | [ ] -> ()
    | x :: tail -> printf "( %7.4f, %7.4f )\n" x.re x.im;
                  complex_show tail


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
    let x = List.map complex_of_int (range 0 (size-1)) in
    let t0 = Sys.time () in
    for i = 1 to size do
        f x;
    done;
    (Sys.time () -. t0) /. (float repeats)


(**************************************************************************************************
 Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.

 Parameters:
   x
     The list of which the DFT will be computed. Given the nature of the implementation, there is no
     restriction on the size of the list, although it will almost always be called with a power of
     two size to give a fair comparison.

 Returns:
     A complex-number list of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let direct_ft x =
    let rec dft_loop x k =                                     (* Compute the kth component; *)
        if k == List.length x then
            [ ]                                                (* FT of an empty vector; *)
        else
            let nx = List.length x in                          (* Length of the vector; *)
            let w = exp { re=0.; im= -2. *. pi *. (float k) /. (float nx) } in
            let wk = List.map (cpow w) (range 0 (nx-1)) in     (* Twiddle factors; *)
            let xwk = List.map2 mul x wk in                    (* Each term of the summation; *)
            let tXk = List.fold_left add zero xwk in           (* Value of the kth component; *)
            tXk :: dft_loop x (k+1)                            (* Recursive call; *)
    in
    (dft_loop x 0)


(**************************************************************************************************
 Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.

 This version is a tail-recursive implementation. It has the same parameters and behaviour of the
 non-tail recursive version, but is probably more time and memory efficient.

 Parameters:
   x
     The list of which the DFT will be computed. Given the nature of the implementation, there is no
     restriction on the size of the list, although it will almost always be called with a power of
     two size to give a fair comparison.

 Returns:
     A complex-number list of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let tr_direct_ft x =
    let rec dft_loop x k tX =
        if k == List.length x then
            tX
        else
            let nx = List.length x in
            let w = exp { re=0.; im= -2. *. pi *. (float k) /. (float nx) } in
            let wk = List.map (cpow w) (range 0 (nx-1)) in
            let xwk = List.map2 mul x wk in
            let tXk = List.fold_left add zero xwk in
            dft_loop x (k+1) (tX @ [ tXk ])
    in
    (dft_loop x 0 [ ])


(**************************************************************************************************
 Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
 complexity.

 Parameters:
   x
     The list of which the DFT will be computed. This should have length power of 2;

 Returns:
   A complex-number list of the same size, with the coefficients of the DFT.
 **************************************************************************************************)
let rec recursive_fft x =
    if List.length x == 1 then
        x
    else
        let rec split x =
            match x with
            | [ ] -> [], []
            | xe :: [ ] -> [ xe ], []
            | xe :: xo :: tail ->
                let xe_tail, xo_tail = split tail in
                ( xe :: xe_tail, xo :: xo_tail )
        in
        let nx = List.length x in
        let w = exp { re=0.; im= -2. *. pi /. (float nx) } in
        let wk = List.map (cpow w) (range 0 ((nx-1)/2)) in
        let (xe, xo) = split x in
        let tXe = recursive_fft xe in
        let tXo = List.map2 mul wk (recursive_fft xo) in
        (List.map2 add tXe tXo) @ (List.map2 sub tXe tXo)


(**************************************************************************************************
 Main function:
 **************************************************************************************************)
let main () =

    (* Start by printing the table with time comparisons: *)
    printf "+---------+---------+---------+---------+---------+---------+\n";
    printf "|    N    |   N^2   | N logN  | Direta  | TRDir.  | Recurs. |\n";
    printf "+---------+---------+---------+---------+---------+---------+\n";

    (* Try it with vectors with size ranging from 32 to 1024 samples: *)
    for r = 5 to 10 do

        (* Compute the average execution time: *)
        let n = int_of_float (2. ** (float r)) in
        let dtime = time_it direct_ft n repeats in
        let tdtime = time_it tr_direct_ft n repeats in
        let rtime = time_it recursive_fft n repeats in

        (* Print the results: *)
        printf "| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f |\n" n (n*n) (n*r) dtime tdtime rtime;
        flush stdout

    done;
    printf "+---------+---------+---------+---------+---------+---------+\n"


let () = main ()
