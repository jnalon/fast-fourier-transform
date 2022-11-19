'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This script compares Cooley-Tukey algorithms.
'' -------------------------------------------------------------------------------------------------

#ifndef __TEST_IT__
#define __TEST_IT__


'' Include necessary libraries:
#include "complex.bas"


'' Definitions:
const Repeats = 500                            '' Number of executions to compute average time;


'' Pretty print of input and output of the Fourier transform for visual inspection.
''
'' Parameters:
''   F
''     Function to be called;
''   Size
''     Number of elements in the vector on which the transform will be applied.
sub TestIt(F as sub(X() as Complex, TX() as Complex), Size as Integer)
    dim X(Size-1) as Complex
    dim TX(Size-1) as Complex
    dim i as Integer
    for i = 0 to Size-1
        X(i) = Complex(i, 0)
    next i
    F(X(), TX())
    print using "N = #### | Input | Output:"; Size
    for i = 0 to Size-1
        print using "  ##: (####.####, ####.####) | (####.####, ####.####)"; _
                    i; X(i).r; X(i).i; TX(i).r; TX(i).i
    next i
    print string(30, "-")
end sub


'' Measure execution time through repeated calls to a (Fast) Fourier Transform function.
''
'' Parameters:
''   f
''     Function to be called;
''   size
''     Size of the vector on which the transform will be applied;
''   repeats
''     Number of times the function will be called. Defaults to REPEAT.
''
'' Returns:
''   The average execution time for that function with a vector of the given size.
function TimeIt(F as sub(X() as Complex, TX() as Complex), Size as Integer, Repeat as Integer = Repeats) as Double
    dim X(Size-1) as Complex
    dim TX(Size-1) as Complex
    dim T0 as Double
    dim i as Integer
    for i = 0 to Size-1
        X(i) = Complex(i, 0)
    next i
    T0 = Timer()
    for i = 1 to Repeat
        F(X(), TX())
    next i
    return (Timer() - T0) / Repeat
end function


#endif
