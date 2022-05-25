/**************************************************************************************************
 * Fast Fourier Transform -- C++ Fixed Point Version
 * This is the header file of the fixed point class and arithmetic.
 **************************************************************************************************/

#ifndef __FIXED_POINT__
#define __FIXED_POINT__


// Include necessary libraries:
#include <iostream>                    // Input and Output;
#include <iomanip>                     // I/O Manipulation;
#include <cmath>                       // Math operations;
#include <vector>                      // Deals with arrays;


using namespace std;


// Definitions:
const int fraction_bits = 16;          ///< Scale factor as a power of two, so it will be faster.
    ///< It's made a constant instead of a property of the fixed point number so numbers with
    ///< different scales don't mix.


/**
 * Define the behavior of fixed point numbers. A number of operators and functions are defined, but
 * they are left without further documentation as what they do is obvious (eg.: operator+ adds two
 * fixed point numbers.
 */
class FixedPoint {
        long int value;                ///< Payload of the fixed point number;
    public:
        FixedPoint() : value(0) { }
        FixedPoint(int value) : value((long int) (value) << fraction_bits) { }
        FixedPoint(long int value) : value(value << fraction_bits) { }
        FixedPoint(float value) : value(value * (1 << fraction_bits)) { }
        FixedPoint(double value) : value(value * (1 << fraction_bits)) { }
        FixedPoint(const FixedPoint &other) : value(other.value) { }
        bool operator <(const FixedPoint &other);
        bool operator <=(const FixedPoint &other);
        FixedPoint operator +(const FixedPoint &other);
        FixedPoint operator -();
        FixedPoint operator -(const FixedPoint &other);
        FixedPoint operator *(const FixedPoint &other);
        FixedPoint operator *(const float &other);
        FixedPoint operator /(const FixedPoint &other);
        FixedPoint operator /(const int &other);
        FixedPoint operator /(const float &other);
        FixedPoint operator %(const FixedPoint &other);
        friend ostream & operator <<(ostream &out, const FixedPoint &number);
};


// Useful constants:
const FixedPoint ZERO = FixedPoint(0);
const FixedPoint PI = FixedPoint(3.14159265359);
const FixedPoint HALFPI = FixedPoint(1.5707963268);
const FixedPoint DOUBLEPI = FixedPoint(6.28318530718);
const vector<FixedPoint> sin_coefficients = {
     2.7557319223985893e-06,   // * x^9 +
     0.0,                      // * x^8 +
    -0.0001984126984126984,    // * x^7 +
     0.0,                      // * x^6 +
     0.008333333333333333,     // * x^5 +
     0.0,                      // * x^4 +
    -0.16666666666666666,      // * x^3 +
     0.0,                      // * x^2 +
     1.0,                      // * x^1 +
     0.0                       // * x^0
};                             ///< Coefficients of the Taylor series for the sine function.


/**
 * Evaluate a polynomial using fixed point numbers. With the correct coefficients coming from a
 * Taylor series, it is possible to compute transcedental functions up to the needed precision. If
 * a polynomial is defined as:
 *
 *   p(x) = a_N x^N + a_{N-1} x^{N-1} + ... + a_2 x^2 + a_1 x + a_0
 *
 * then the coefficients vector should be of the form:
 *
 *   C = [ a_N, a_{N-1}, ..., a_2, a_1, a_0 ]
 *
 * @param x `FixedPoint` number of which the polynomial will be computed.
 * @param coefficients vector of `FixedPoint` coefficients as described above.
 * @return The result of the computation.
 */
FixedPoint evaluate_polynomial(FixedPoint x, vector<FixedPoint> coefficients);


/**
 * Compute the sine of a number given in radians and in the `FixedPoint` format.
 *
 * @param x `FixedPoint` number in radians of which the sine will be computed.
 * @return The result of the operation.
 */
FixedPoint sin(FixedPoint x);


/**
 * Compute the cosine of a number given in radians and in the `FixedPoint` format.
 *
 * @param x `FixedPoint` number in radians of which the cosine will be computed.
 * @return The result of the operation.
 */
FixedPoint cos(FixedPoint x);


#endif
