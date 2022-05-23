/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is a header file for a simple library of complex numbers.
 **************************************************************************************************/

#ifndef __MY_COMPLEX__
#define __MY_COMPLEX__


// Include necessary libraries:
#include <stdlib.h>                            // Standard Library;
#include <stdio.h>                             // Input and Output;
#include <math.h>                              // Math Functions;


/**
 Data structure to deal with complex numbers. There is, of course, an external library that
 implements operations for comple numbers, but I'm defining this here to reduce dependencies.
 */
typedef struct {
    float r;                                   //< Real part of the complex number;
    float i;                                   //< Imaginary part of the complex number;
} Complex;


Complex cmplx(float r, float i);               //< Initialization of a complex number.
Complex cadd(Complex a, Complex b);            //< Complex addition of numbers a and b.
Complex csub(Complex a, Complex b);            //< Complex subtraction of number b from a.
Complex cmul(Complex a, Complex b);            //< Complex product of numbers a and b.
Complex cexpn(float angle);                    //< Complex exponential of an angle.


#endif
