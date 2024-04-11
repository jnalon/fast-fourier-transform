/**************************************************************************************************
 * Fast Fourier Transform -- C# Version
 * This is the implementation of the complex number library
 **************************************************************************************************/

// Include necessary libraries:
using System;                                  // Input and output and standard library;


public class Complex
{
    public double re;                          // Real part;
    public double im;                          // Imaginary part;

    public Complex()
    {
        this.re = 0;
        this.im = 0;
    }

    public Complex(double real, double imag)
    {
        this.re = real;
        this.im = imag;
    }

    // Sum of two complex numbers:
    public static Complex operator +(Complex z1, Complex z2)
    {
        return new Complex(z1.re + z2.re, z1.im + z2.im);
    }

    // Difference of two complex numbers:
    public static Complex operator -(Complex z1, Complex z2)
    {
        return new Complex(z1.re - z2.re, z1.im - z2.im);
    }

    // Product of two complex numbers:
    public static Complex operator *(Complex z1, Complex z2)
    {
        return new Complex(z1.re*z2.re - z1.im*z2.im, z1.re*z2.im + z1.im*z2.re);
    }

    // Product with a scalar:
    public static Complex operator *(Complex z, double x)
    {
        return new Complex(z.re*x, z.im*x);
    }

    // Complex Exponential:
    public static Complex exp(double a)
    {
        return new Complex(Math.Cos(a), Math.Sin(a));
    }

    // String representation for printing:
    public override string ToString()
    {
        return(String.Format("{0} + {1}i", re, im));
    }
}
