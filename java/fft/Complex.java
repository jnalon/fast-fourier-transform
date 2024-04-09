/**************************************************************************************************
 * Fast Fourier Transform -- Java Version
 * This is the implementation of the complex number library.
 **************************************************************************************************/

package fft;


public class Complex {

    public float r, i;

    // Constructor:
    public Complex(float re, float im)
    {
        r = re;
        i = im;
    }

    // Constructor:
    public Complex()
    {
        this(0, 0);
    }

    // Add the argument to this, giving the result as a new complex number:
    public Complex add(Complex c)
    {
        return new Complex(r + c.r, i + c.i);
    }

    // Subtract the argument from this, giving the result as a new complex number:
    public Complex sub(Complex c)
    {
        return new Complex(r - c.r, i - c.i);
    }

    // Multiply the argument with this, giving the result as a new complex number:
    public Complex mul(Complex c)
    {
        return new Complex(r*c.r - i*c.i, r*c.i + i*c.r);
    }

    // Divide this by the argument, giving the result as a new complex number:
    public Complex div(float a)
    {
        return new Complex(r/a, i/a);
    }

    // Complex exponential of an angle:
    public static Complex exp(float a)
    {
        return new Complex((float) Math.cos(a), (float) Math.sin(a));
    }

}
