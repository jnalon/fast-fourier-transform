/**************************************************************************************************
 * Fast Fourier Transform -- C++ Fixed Point Version
 * This is the implementation file of the fixed point class and arithmetic.
 **************************************************************************************************/

// Include necessary libraries:
#include "fixed_point.h"


bool FixedPoint::operator<(const FixedPoint &other) {
    return value < other.value;
}

bool FixedPoint::operator<=(const FixedPoint &other) {
    return value <= other.value;
}

FixedPoint FixedPoint::operator+(const FixedPoint &other) {
    FixedPoint c;
    c.value = value + other.value;
    return c;
}

FixedPoint FixedPoint::operator-() {
    FixedPoint c;
    c.value = - value;
    return c;
}

FixedPoint FixedPoint::operator-(const FixedPoint &other) {
    FixedPoint c;
    c.value = value - other.value;
    return c;
}

// The function below works great and fast because we have access to long ints. If they weren't
// available, we would need to do something like (default_scale = 1 << fraction_bits):
//     int i1 = value / default_scale;
//     int f1 = value % default_scale;
//     int i2 = other.value / default_scale;
//     int f2 = other.value % default_scale;
//     int i3 = (i1 * i2) * default_scale;
//     int f3 = (i1 * f2 + i2 * f1) + (f1 * f2) / default_scale;
//     return FixedPoint(i3 + f3);
FixedPoint FixedPoint::operator*(const FixedPoint &other) {
    FixedPoint c;
    c.value = (value * other.value) >> fraction_bits;
    return c;
}

FixedPoint FixedPoint::operator*(const float &other) {
    return (*this) * FixedPoint(other);
}

FixedPoint FixedPoint::operator/(const FixedPoint &other) {
    FixedPoint c;
    c.value = (value << fraction_bits) / other.value;
    return c;
}

FixedPoint FixedPoint::operator/(const int &other) {
    FixedPoint c;
    c.value = value / other;
    return c;
}

FixedPoint FixedPoint::operator/(const float &other) {
    return (*this) / FixedPoint(other);
}

FixedPoint FixedPoint::operator%(const FixedPoint &other) {
    FixedPoint c;
    c.value = value % other.value;
    return c;
}

ostream &operator<<(ostream &out, const FixedPoint &number) {
    out << ((float) number.value / (float) (1 << fraction_bits));
    return out;
}


FixedPoint evaluate_polynomial(FixedPoint x, vector<FixedPoint> coefficients) {
    FixedPoint total(0);
    for(FixedPoint c : coefficients)
        total = x * total + c;
    return total;
}


FixedPoint sin(FixedPoint x) {
    if (x < ZERO)
        return - sin(-x);
    if (x <= HALFPI)
        return evaluate_polynomial(x, sin_coefficients);
    else if (x <= PI) {
        FixedPoint argument = - x + PI;
        return evaluate_polynomial(argument, sin_coefficients);
    } else {
        FixedPoint argument = (x + PI) % DOUBLEPI - PI;
        return sin(argument);
    }
}


FixedPoint cos(FixedPoint x) {
    return sin(x + HALFPI);
}
