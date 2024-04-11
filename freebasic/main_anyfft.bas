'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This script compares Cooley-Tukey algorithms for powers of 2 only.
'' -------------------------------------------------------------------------------------------------

'' Import necessary libraries:
#include "test_it.bas"
#include "fft.bas"


dim SIZES(1 to 8) as Integer => { 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 }
dim n as Integer, i as Integer
dim double_directft_time as Double, double_recursiveft_time as Double

print "+---------+---------+---------+---------+"
print "|    N    |   N^2   | Direct  | Recurs. |"
print "+---------+---------+---------+---------+"

'' Try it with vectors with size ranging from 32 to 1024 samples:
for i = 1 to 8

    '' Compute the average execution time:
    n = SIZES(i)
    double_directft_time = TimeIt(@DirectFT, n)
    double_recursiveft_time = TimeIt(@RecursiveNFFT, n)

    '' Print the results:
    print using "| ####### | ####### | ##.#### | ##.#### |"; _
                n; n*n; _
                double_directft_time; double_recursiveft_time

next i

print "+---------+---------+---------+---------+"
