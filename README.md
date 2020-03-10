# Fast Fourier Transform
## Simple implementations of the Fast Fourier Transform in various languages

This project was created as a demonstration of my personal effort to learn programming languages. When I'm learning a new language, I try a certain number of algorithms to see how much I was able to understand, and how much time it takes until I can produce code. Among the algorithms I always try are: Fibonacci sequence, Factorial, Quicksort, The Knight's Round, Tower of Hanoi, Fast Fourier Transform, and others. Since I'm an Electronic Engineer with heavy training in Computers and Telecomunications, I chose the FFT to investigate with more details.

Every implementation here does basically the same: they execute the Cooley-Tukey decomposition FFT for vectors of various lengths, almost always a power of two, although some implementations can deal with composite numbers. The functions are executed a number of times, and execution time is measured and averaged. A table with a comparison of efficiency is printed on screen.

Currently, this repository contains implementations in the following languages:

* To be added

I don't claim that any of these implementations is the most efficient way to deal with the FFT in any language. In fact, most likely, it isn't, as most of them have bindings for a tested and fast library, such as FFTW.
