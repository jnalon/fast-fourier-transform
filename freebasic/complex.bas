'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This file implements a very simple library of complex numbers.
'' -------------------------------------------------------------------------------------------------

#ifndef __COMPLEX__
#define __COMPLEX__


'' Definition of PI, since it's not pre-defined:
const PI as Double = 4 * atn(1)


'' Data structure to deal with complex numbers.
type Complex
    as Double r                                '' Real part of the complex number;
    as Double i                                '' Imaginary part of the complex number;
    declare constructor
    declare constructor(as Double, as Double)
    declare operator cast() as String
    declare static function exp(as Double) as Complex
end type

'' Initializes a complex number with zeros:
constructor Complex
    r = 0.0
    i = 0.0
end constructor

'' Initializes a complex number with a real and imaginary part:
constructor Complex(real as Double, imag as Double=0.0)
    r = real
    i = imag
end constructor

'' String representation of the number:
operator Complex.cast() as String
    return "(" & r & ", " & i & ")"
end operator

'' Complex exponential of an angle:
function Complex.exp(angle as Double) as Complex
    return Complex(cos(angle), sin(angle))
end function

'' Complex addition:
operator + (c1 as Complex, c2 as Complex) as Complex
    return Complex(c1.r + c2.r, c1.i + c2.i)
end operator

'' Complex subtraction:
operator - (c1 as Complex, c2 as Complex) as Complex
    return Complex(c1.r - c2.r, c1.i - c2.i)
end operator

'' Complex multiplication:
operator * (c1 as Complex, c2 as Complex) as Complex
    return Complex(c1.r*c2.r - c1.i*c2.i, c1.r*c2.i + c1.i*c2.r)
end operator

'' Complex multiplication by a scalar:
operator * (a as Double, c as Complex) as Complex
    return Complex(a * c.r, a * c.i)
end operator

'' Division of a Complex number by a scalar:
operator / (c as Complex, a as Double) as Complex
    return Complex(c.r/a, c.r/a)
end operator


'' Auxiliary method to print an array of complex numbers:
sub ComplexShow(x() as Complex, prefix as String = "")
    dim N as Integer = UBound(x)
    dim i as Integer
    for i = 0 to N
        Print prefix; X(i)
    next i
end sub


#endif
