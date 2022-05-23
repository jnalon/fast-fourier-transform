/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the implementation file for the complex number library.
 **************************************************************************************************/

// Include necessary libraries:
#include "my_complex.h"


Complex cmplx(float r, float i)
{
    Complex w;

    w.r = r;                                   // Real part;
    w.i = i;                                   // Imaginary part;
    return w;
}


Complex cadd(Complex a, Complex b)
{
    Complex w;

    w.r = a.r + b.r;                           // Real part of sum;
    w.i = a.i + b.i;                           // Imaginary part of sum;
    return w;
}


Complex csub(Complex a, Complex b)
{
    Complex w;

    w.r = a.r - b.r;                           // Real part of difference;
    w.i = a.i - b.i;                           // Imaginary part of difference;
    return w;
}


Complex cmul(Complex a, Complex b)
{
    Complex w;

    w.r = a.r*b.r - a.i*b.i;                   // Real part of product;
    w.i = a.r*b.i + a.i*b.r;                   // Imaginary part of product;
    return w;
}


Complex cexpn(float angle)
{
    Complex w;

    w.r = cos(angle);                          // Real part of exponential;
    w.i = sin(angle);                          // Imaginary part of exponential;
    return w;
}
