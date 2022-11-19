'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This script compares Cooley-Tukey algorithms for powers of 2 only.
'' -------------------------------------------------------------------------------------------------

'' Import necessary libraries:
#include "test_it.bas"
#include "fft.bas"


dim r as Integer, n as Integer
dim double_directft_time as Double, double_recursiveft_time as Double, _
    double_iterativeft_time as Double

print "+---------+---------+---------+---------+---------+---------+"
print "|    N    |   N^2   | N log N | Direct  | Recurs. | Itera.  |"
print "+---------+---------+---------+---------+---------+---------+"

'' Try it with vectors with size ranging from 32 to 1024 samples:
for r = 5 to 10

    '' Compute the average execution time:
    n = 2 ^ r
    double_directft_time = TimeIt(@DirectFT, n)
    double_recursiveft_time = TimeIt(@RecursiveFFT, n)
    double_iterativeft_time = TimeIt(@IterativeFFT, n)

    '' Print the results:
    print using "| ####### | ####### | ####### | ##.#### | ##.#### | ##.#### |"; _
                n; n*n; r*n; _
                double_directft_time; double_recursiveft_time; double_iterativeft_time

next r

print "+---------+---------+---------+---------+---------+---------+"
