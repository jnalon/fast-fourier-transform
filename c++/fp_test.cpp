/**************************************************************************************************
 * Fast Fourier Transform -- C++ Fixed Point Version
 * Test program for FixedPoint and Complex libraries.
 **************************************************************************************************
 * This program doesn't need much to be compiled and run. It can be done, as far as I know, with
 * any C++ compiler, just remember to link the math library. In my box, I used the command:
 *
 * $ g++ -o fp_test fp_test.cpp fixed_point.cpp -lm
 *
 * It can be run with the command (remember to change permission to execute):
 *
 * $ ./fp_test
 **************************************************************************************************/

// Include necessary libraries:
#include <iostream>                    // Input and Output;
#include <iomanip>                     // I/O Manipulation;
#include <cmath>                       // Math operations;
#include <vector>                      // Deals with arrays;
#include "fixed_point.h"               // Fixed Point library;
#include "my_complex.h"                // Complex numbers;


/**************************************************************************************************
 Main Function:
 **************************************************************************************************/
int main(int argc, char *argv[]) {

    FixedPoint a(3.141592), b(1.570796), c(2);
    std::cout << "A = " << a << std::endl;
    std::cout << "B = " << b << std::endl;
    std::cout << "C = " << c << std::endl;
    std::cout << "A + B = " << a + b << std::endl;
    std::cout << "A - B = " << a - b << std::endl;
    std::cout << "B - A = " << b - a << std::endl;
    std::cout << "A * B = " << a * b << std::endl;
    std::cout << "A / B = " << a / b << std::endl;
    std::cout << "B / A = " << b / a << std::endl;
    std::cout << "A / C = " << a / c << std::endl;
    std::cout << "B / C = " << b / c << std::endl;
    std::cout << "A / 2 = " << a / 2.0f << std::endl;
    std::cout << "A / 2 = " << a / 2 << std::endl;

    vector<FixedPoint> coefficients = { 1.0, 2.0, 3.0 };
    std::cout << "C = " << c << std::endl;
    std::cout << "Coefficients = ";
    for(auto c : coefficients)
        std::cout << c << " ";
    std::cout << std::endl;
    std::cout << "C^2 + 2*C + 3 = " << evaluate_polynomial(c, coefficients);
    std::cout << std::endl;

    std::cout << "sin(A) = " << sin(a) << std::endl;
    std::cout << "sin(B) = " << sin(b) << std::endl;

    for(int i=-16; i<=16; i++) {
        FixedPoint x = i * 3.14159265359 / 8;
        std::cout << x << " = ";
        std::cout << (x + PI) % DOUBLEPI - PI << " = ";
        std::cout << sin(x) << " | " << cos(x) << std::endl;
    }

    Complex<float> z;
    Complex<float> w(1.0f, 2.0f);
    Complex<float> v(2.1f, 1.2f);

    std::cout << "z = " << z << std::endl;
    std::cout << "w = " << w << std::endl;
    std::cout << "v = " << v << std::endl;

    std::cout << "z + w = " << z + w << std::endl;
    std::cout << "w + v = " << w + v << std::endl;
    std::cout << "w - v = " << w - v << std::endl;
    std::cout << "v - w = " << v - w << std::endl;
    std::cout << "w * v = " << w * v << std::endl;

    std::cout << "cexp(2) = " << cexpn(2.0f) << std::endl;
    std::cout << "cexp(c) = " << cexpn(c) << std::endl;

    vector<Complex<float>> x(16);
    for(int i=0; i<16; i++)
        x[i] = i + 1;

    for(int i=0; i<16; i++)
        std::cout << x[i] << std::endl;
}
