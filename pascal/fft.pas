program FFT;

{--------------------------------------------------------------------------------------------------
 Fast Fourier Transform -- Pascal Version
 This version implements Cooley-Tukey algorithm for powers of 2 only.

 José Alexandre Nalon
 --------------------------------------------------------------------------------------------------
 This version of the algorithm was compiled and tested using FreePascal, but it should work with
 other compilers. You can compile using the following command:

 $ fpc fft.pas

 And run by issuing (remember to give permission to execute):

 $ ./fft
 --------------------------------------------------------------------------------------------------}

{ Using needed libraries. }
uses
    dos,                               { Timing functions; }
    math;                              { Math functions; }


{--------------------------------------------------------------------------------------------------
 Definitions
 --------------------------------------------------------------------------------------------------}
{ Complex number arithmetics. Many of the Pascal compilers out there have a library for complex
  numbers, but they don't behave the same, so I made my own.                   }
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
   Measure execution time through repeated calls to a (Fast) Fourier Transform function.

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
    x: ComplexArray;
    h, m, s, c: word;
    t0, t1: real;
    j: integer;
begin
    for j := 0 to size-1 do                    { Initialize the vector and transform; }
        x[j] := j;
    GetTime(h, m, s, c);                       { Start counting time; }
    t0 := h*3600 + m*60 + s + c/100;
    for j := 1 to repeats do
        f(x, size);
    GetTime(h, m, s, c);                       { Finish counting time; }
    t1 := h*3600 + m*60 + s + c/100;
    TimeIt := (t1 - t0) / repeats;
end;


{--------------------------------------------------------------------------------------------------
 Function: DirectFT
   Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
   complexity.

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
    tX: ComplexArray;                          { Temporary results; }
    W, Wk, Wkn: Complex;                       { Twiddle factors; }
    j, k: integer;
begin
    W := cexp(-2*Pi/n);                        { Initialize twiddle factors; }
    Wk := Complex(1);
    for k := 0 to n-1 do begin
        tX[k] := 0;                            { Accumulate the results; }
        Wkn := Complex(1);                     { Initialize twiddle factors; }
        for j := 0 to n-1 do begin
            tX[k] := tX[k] + Wkn * x[j];
            Wkn := Wkn * Wk;                   { Update twiddle factors; }
        end;
        Wk := Wk * W;
    end;
    DirectFT := tX;                            { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Function: RecursiveFFT
   Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
   complexity.

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
    xe, xo, tXe, tXo, tX: ComplexArray;        { Vectors with intermediate results; }
    W, Wk, wt: Complex;                        { Twiddle factors; }
    k, n2: integer;
begin
    if n=1 then                                { A length-1 vector is its own FT; }
        tX[0] := x[0]
    else begin
        n2 := n div 2;

        for k := 0 to n2-1 do begin            { Split even and odd samples; }
            xe[k] := x[2*k];
            xo[k] := x[2*k+1];
        end;
        tXe := RecursiveFFT(xe, n2);           { Transform of even samples; }
        tXo := RecursiveFFT(xo, n2);           { Transform of odd samples; }

        W := cexp(-2*Pi/n);
        Wk := Complex(1);
        for k := 0 to n2-1 do begin
            wt := Wk * tXo[k];                 { Recombine results; }
            tX[k] := tXe[k] + wt;
            tX[k+n2] := tXe[k] - wt;
            Wk := Wk * W;                      { Update twiddle factors; }
        end;

    end;
    RecursiveFFT := tX;                        { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Function: BitReverse
   Bit-reversed version of an integer number.

 Parameters:
   k
     The number to be bit-reversed;
   r
     The number of bits to take into consideration when reversing.

  Returns:
   The number k, bit-reversed according to integers with r bits.
 --------------------------------------------------------------------------------------------------}
function BitReverse(k, r: integer): integer;
var
    l, i: integer;
begin
    l := 0;                                    { Accumulate the results; }
    for i := 1 to r do begin                   { Loop on every bit; }
        l := l * 2 + (k mod 2);
        k := k div 2;                          { Test next bit; }
    end;
    BitReverse := l;                           { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Function: IterativeFFT
   Fast Fourier Transform using an iterative in-place decimation in time algorithm. This has
   O(N log_2(N)) complexity, and since there are less function calls, it will probably be
   marginally faster than the recursive versions.

 Parameters:
   x
     The vector of which the FFT will be computed. This should always be called with a vector of
     a power of two length, or it will fail. No checks on this are made.
   n
     The number of elements in the vector.

 Returns:
    A complex-number vector of the same size, with the coefficients of the DFT.
 --------------------------------------------------------------------------------------------------}
function IterativeFFT(var x: ComplexArray; n: integer): ComplexArray;
var
    tX: ComplexArray;                          { Temporary results; }
    W, Wkn: Complex;                           { Twiddle factors; }
    r, k, l, m, p, q, step: integer;
begin
    r := Round(log2(n));                       { Number of bits; }
    for k := 0 to n-1 do begin
        l := BitReverse(k, r);                 { Reorder the vector according to the }
        tX[l] := x[k];                         {   bit-reversed order; }
    end;

    step := 1;                                 { Auxiliary for computation of twiddle factors; }
    for k := 0 to r-1 do begin
        l := 0;
        while l < n-1 do begin
            W := cexp(-Pi/step);               { Twiddle factors; }
            Wkn := Complex(1);
            for m := 0 to step-1 do begin
                p := l + m;
                q := p + step;
                tX[q] := tX[p] - Wkn * tX[q];  { Recombine results; }
                tX[p] := 2*tX[p] - tX[q];
                Wkn := Wkn * W;                { Update twiddle factors; }
            end;
            l := l + 2*step;
        end;
        step := 2*step;
    end;

    IterativeFFT := tX;                         { Return value; }
end;


{--------------------------------------------------------------------------------------------------
 Main Program:
 --------------------------------------------------------------------------------------------------}
var
    f: DFTFunction;
    dtime, rtime, itime: real;
    r, n: integer;
begin
    // Start by printing the table with time comparisons:
    writeln('+---------+---------+---------+---------+---------+---------+');
    writeln('|    N    |   N^2   | N logN  | Direct  | Recurs. | Inter.  |');
    writeln('+---------+---------+---------+---------+---------+---------+');

    // Try it with vectors with size ranging from 32 to 1024 samples:
    for r := 5 to 10 do begin

        // Compute the average execution time:
        n := floor(power(2, r));
        f := @DirectFT;     dtime := timeit(f, n, Repeats);
        f := @RecursiveFFT; rtime := timeit(f, n, Repeats);
        f := @IterativeFFT; itime := timeit(f, n, Repeats);

        // Print the results:
        writeln('| ', n:7, ' | ', n*n:7, ' | ', r*n:7, ' | ', dtime:7:4, ' | ', rtime:7:4, ' | ', itime:7:4, ' |');
    end;
    writeln('+---------+---------+---------+---------+---------+---------+');
end.
