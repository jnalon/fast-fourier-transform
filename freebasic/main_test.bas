'' -------------------------------------------------------------------------------------------------
'' Fast Fourier Transform - FreeBasic Vresion
'' This script compares Cooley-Tukey algorithms.
'' -------------------------------------------------------------------------------------------------

'' Include necessary libraries:
#include "fft.bas"
#include "test_it.bas"


print "Direct FT -"
TestIt(@DirectFT, 8)
print "Recursive FFT -"
TestIt(@RecursiveFFT, 8)
print "Iterative FFT -"
TestIt(@IterativeFFT, 8)
print "Direct FT -"
TestIt(@DirectFT, 16)
print "Recursive FFT -"
TestIt(@RecursiveFFT, 16)
print "Iterative FFT -"
TestIt(@IterativeFFT, 16)

print "Direct FT -"
TestIt(@DirectFT, 12)
print "Recursive FFT -"
TestIt(@RecursiveNFFT, 12)
