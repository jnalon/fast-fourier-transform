/**************************************************************************************************
 * Fast Fourier Transform -- C++ Version
 * This is the header file of the complex class and arithmetic.
 **************************************************************************************************/

#ifndef __MY_COMPLEX__
#define __MY_COMPLEX__


// Include necessary libraries:
#include <iostream>                    // Input and Output;
#include <iomanip>                     // I/O Manipulation;
#include <cmath>                       // Math operations;


using namespace std;


/**
 * Small class to operate with complex numbers. A small number of operations are defined, but they
 * are not documented as they behavior is obvious. Eg.: operator+ adds two complex numbers. Real and
 * imaginary parts are published for easy access.
 */
template<typename T>
class Complex {
    public:
        T r;                  //< Real part of the complex number;
        T i;                  //< Imaginary part of the complex number;
        Complex<T>() : r(0), i(0) { }
        Complex<T>(T re, T im=0) : r(re), i(im) { }
        Complex<T> operator +(Complex<T> c) {
            return Complex(r + c.r, i + c.i);
        }
        Complex<T> operator -(Complex<T> c) {
            return Complex(r - c.r, i - c.i);
        }
        Complex<T> operator *(Complex<T> c) {
            return Complex(r * c.r - i * c.i, r * c.i + i * c.r);
        }
        Complex<T> operator *(T a) {
            return Complex(r * a, i * a);
        }
};


/**
 * Overload of the stream << operator to deal with Complex numbers of any type.
 */
template<typename T>
ostream & operator <<(ostream &out, const Complex<T> &z) {
    out << "(" << z.r << ", " << z.i << ")";
    return out;
}


/**
 * Compute the complex exponential of an angle given in radians. The result is a complex number
 * given by
 *
 *   e^(ja) = cos(a) + j sin(a)
 *
 * @param a angle in radians to be computed.
 * @return The result of the operation as described above.
 */
template<typename T>
Complex<T> cexpn(T a) {
    return Complex<T>(cos(a), sin(a));
}


#endif
