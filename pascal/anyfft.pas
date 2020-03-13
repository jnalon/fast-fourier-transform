program FFT;

{--------------------------------------------------------------------------------------------------
 Fast Fourier Transform -- Pascal Version
 This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).

 José Alexandre Nalon
 --------------------------------------------------------------------------------------------------
 This version of the algorithm was compiled and tested using FreePascal, but it should work with
 other compilers. You can compile using the following command:

 $ fpc fft.pas

 And run by issuing (remember to give permission to execute):

 $
 --------------------------------------------------------------------------------------------------}

{ Using needed libraries. }
uses
    dos,                               { Timing functions; }
    math;                              { Math functions; }


{--------------------------------------------------------------------------------------------------
 Definitions
 --------------------------------------------------------------------------------------------------}
{ Class that implements complex number arithmetics. Many of the Pascal compilers out there have a
  library for complex numbers, but they don't behave the same, so I made my own.                   }
type
    Complex = record
        r: real;
        i: real;
    end;

{ Vectors of complex numbers, for simplicity, assume they will be 1024 samples or less. }
type ComplexArray = array [0..1023] of Complex;

{ Prototype for a DFT computation procedure.                                                       }
type DFTFunction = function(var x: ComplexArray; n: integer): ComplexArray;


{--------------------------------------------------------------------------------------------------
 Definitions:
 --------------------------------------------------------------------------------------------------}
const
    Repeats = 500;                     { Number of executions to compute average time; }


{--------------------------------------------------------------------------------------------------
 Small set of functions to operate with complex numbers.
 --------------------------------------------------------------------------------------------------}
operator := (r: real) z: Complex;                      { Casting from real numbers; }
begin
    z.r := r;
    z.i := 0.;
end;

operator + (a, b: Complex) z: Complex;                 { Complex addition of numbers a and b; }
begin
    z.r := a.r + b.r;
    z.i := a.i + b.i;
end;

operator - (a, b: Complex) z: Complex;                 { Complex subtraction of number b from a; }
begin
    z.r := a.r - b.r;
    z.i := a.i - b.i;
end;

operator * (a, b: Complex) z: Complex;                 { Complex product of numbers a and b; }
begin
    z.r := a.r*b.r - a.i*b.i;
    z.i := a.r*b.i + a.i*b.r;
end;

operator * (a: Real; b: Complex) z: Complex;           { Product with an scalar; }
begin
    z.r := a * b.r;
    z.i := a * b.i;
end;

function cexp(r: Real): Complex;                       { Complex exponential of an angle; }
begin
    cexp.r := cos(r);
    cexp.i := sin(r);
end;


{--------------------------------------------------------------------------------------------------
 Auxiliary procedure: ComplexShow
   Pretty printing of an array of complex numbers, used to inspect results.

 Parameters:
   x
     A vector of complex numbers, according to the definition above;
   n
     Number of elements on the vector.
 --------------------------------------------------------------------------------------------------}
procedure ComplexShow(var x: ComplexArray; n: integer);
var
    i: integer;
begin
    for i := 0 to n-1 do
        WriteLn('(', x[i].r:7:4, ', ', x[i].i:7:4, ')');
end;


{--------------------------------------------------------------------------------------------------
 Auxiliary function: TimeIt
   This function calls a Fast Fourier Transform function repeatedly a certain number of times,
   measure execution time and average it.

 Parameters:
  f
    Function to be called, with the given prototype. The first complex vector is the input
    vector, the second complex vector is the result of the computation, and the integer is the
    number of elements in the vector;
  size
    Number of elements in the vector on which the transform will be applied;
  repeat
    Number of times the function will be called.

 Returns:
   The average execution time for that function with a vector of the given size.
 --------------------------------------------------------------------------------------------------}
function TimeIt(var f: DFTFunction; size, repeats: integer): real;
var
    x, y: ComplexArray;
    h, m, s, c: word;
    t0, t1: real;
    j: integer;
begin
    for j := 0 to size-1 do                    { Initializes the vector and transform; }
        x[j] := j;
    GetTime(h, m, s, c);                       { Starts counting time; }
    t0 := h*3600 + m*60 + s + c/100;
    for j := 1 to repeats do
        y := f(x, size);
    GetTime(h, m, s, c);                       { Finishes counting time; }
    t1 := h*3600 + m*60 + s + c/100;
    TimeIt := (t1 - t0) / repeats;
end;


{--------------------------------------------------------------------------------------------------
 Function: DirectFT
   Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that has
   O(N^2) complexity.

 Parameters:
   x
     The vector of which the DFT will be computed. Given the nature of the implementation, there
     is no restriction on the size of the vector, although it will almost always be called with a
     power of two size to give a fair comparison;
   n
     The number of elements in the vector.

 Returns:
    A complex-number vector of the same size, with the coefficients of the DFT.
 --------------------------------------------------------------------------------------------------}
function DirectFT(var x: ComplexArray; n: integer): ComplexArray;
var
    y: ComplexArray;                           { Temporary results; }
    W, Wk, Wkn: Complex;                       { Twiddle factors; }
    j, k: integer;
begin
    W := cexp(-2*Pi/n);                        { Initializes twiddle factors; }
    Wk := Complex(1);
    for k := 0 to n-1 do begin
        y[k] := 0;                             { Accumulates the results; }
        Wkn := Complex(1);                     { Initializes twiddle factors; }
        for j := 0 to n-1 do begin
            y[k] := y[k] + Wkn * x[j];
            Wkn := Wkn * Wk;                   { Updates twiddle factors; }
        end;
        Wk := Wk * W;
    end;
    DirectFT := y;                             { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Function: Factor
   This function finds the smallest prime factor of a given number. If the argument is prime
   itself, then it is the return value.

 Parameters:
   n
     Number to be inspected.

 Returns:
   The smallest prime factor, or the number itself if it is already a prime.
 --------------------------------------------------------------------------------------------------}
function factor(n: integer): integer;
var
    rn, i: integer;
begin
    rn := floor(sqrt(n));                      { Search up to the square root of the number; }
    factor := n;
    for i := 2 to rn do
        if (n mod i) = 0 then begin
            factor := i;                       { If remainder is zero, a factor is found; }
            break;
        end;
end;


{--------------------------------------------------------------------------------------------------
 Function: RecursiveFFT
   Computes the Fast Fourier Ttransform using a recursive decimation in time algorithm. This has
   smallest complexity than the direct FT, though the exact value is difficult to compute.

 Parameters:
   x
     The vector of which the FFT will be computed. This should always be called with a vector of
     a power of two length, or it will fail. No checks on this are made.
   n
     The number of elements in the vector.

 Returns:
    A complex-number vector of the same size, with the coefficients of the DFT.
 --------------------------------------------------------------------------------------------------}
function RecursiveFFT(var x: ComplexArray; n: integer): ComplexArray;
var
    xj, FXj, y: ComplexArray;                  { Vectors with intermediate results; }
    W, Wj, Wkj: Complex;                       { Twiddle factors; }
    j, k, n1, n2: integer;
begin
    n1 := factor(n);                           { Smallest prime factor of length; }
    if n1=n then                               { If the length is prime itself, }
        y := DirectFT(x, n)                    {   the transform is given by the direct form; }
    else begin
        for k := 0 to n-1 do                   { Initializes results with 0; }
            y[k] := 0;
        n2 := Floor(n / n1);                   { Decompose in two factors, N1 being prime; }
        W := cexp(-2*Pi/n);                    { Twiddle factor;}
        Wj := Complex(1);
        for j := 0 to n1-1 do begin            { Computes every subsequence of size N2; }
            for k := 0 to N2-1 do
                xj[k] := x[k*n1+j];            { Creates the subsequence; }
            FXj := RecursiveFFT(xj, n2);       { Computes the DFT of the subsequence; }
            Wkj := Complex(1);
            for k := 0 to N-1 do begin
                y[k] := y[k] + FXj[k mod n2] * Wkj;
                Wkj := Wkj * Wj;               { Update twiddle factors; }
            end;
            Wj := Wj * W;
        end;
    end;
    RecursiveFFT := y;                         { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Main Program:
 --------------------------------------------------------------------------------------------------}
var
    f: DFTFunction;
    sizes: array[0..7] of integer = (2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5);
    dtime, rtime: real;
    n, i: integer;
begin
    // Starts by printing the table with time comparisons:
    writeln('+---------+---------+---------+---------+');
    writeln('|    N    |   N^2   | Direta  | Recurs. |');
    writeln('+---------+---------+---------+---------+');

    // Try it with vectors with the given sizes:
    for i := 0 to 7 do begin

        // Computes the average execution time:
        n := sizes[i];
        f := @DirectFT;     dtime := timeit(f, n, Repeats);
        f := @RecursiveFFT; rtime := timeit(f, n, Repeats);

        // Print the results:
        writeln('| ', n:7, ' | ', n*n:7, ' | ', dtime:7:4, ' | ', rtime:7:4, ' |');
    end;
    writeln('+---------+---------+---------+---------+');
end.
