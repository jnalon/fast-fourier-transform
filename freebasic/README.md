# FreeBASIC Version
This folder contains the FreeBASIC version of the Discrete Fourier Transform. If, like me, you were a kid playing with and learning about computers in the 80's, you are very familiar with the BASIC language. Every home computer of that time had its own version of the language, and we marveled at how we could issue simple commands and the computer would do what we told it. Most people of that time had their start with BASIC, until computers got a little more powerful and other languages were made available.

BASIC (Beginner's All-purpose Symbolic Instruction Code) was created by John G. Kemeny and Thomas E. Kurtz in 1964, as a programming language aimed at non-technical people, with simplified ideas from FORTRAN and Algol. It was an awful language, where all variables were global, without structured constructs and with very limited subroutines and functions support. Probably exactly because of its simplicity, it became the norm in home computers in the 80's.

Its popularity declined in the beginning of the 90's, as more powerful computers were created, which were able to give support to more modern languages. The language, however, was never actually abandoned, and evolved to become better. Microsoft released QuickBASIC, which implemented a lot of structured constructs such as while loops, user-defined types, and actual subroutines. A lot of these extensions were inspired by C++ and Pascal, which were popular languages at that time. FreeBASIC is a free open source implementation compatible with QuickBASIC, but expanding on it. You can find more about it on their [website](https://www.freebasic.net/).

## Comments on the Language
FreeBASIC allows for user-defined types, functions, subroutines and structured loops. As such, it isn't very different from Pascal or similar languages. User types can be used to program classes and objects, and it is very easy to split a complex program in modules for easy maintenance.

Operator overloading is available, which allows to simplify expressions. However, instead of making the operator a method of a class, it is must be defined outside its scope, taking the class as the types of its arguments. I imagine this makes implementation and compilation easier, and it mimics a possible behaviour from C++. While this is not a problem, code would probably be a little cleaner if it was the other way around.

Returning arrays from functions is not easy, so if you need to compute an array, its descriptor must be passed to a subroutine as an argument. This means that all computations are made in place, which can cause side effects if care is not taken.

Also, FreeBASIC is not case sensitive, which means that, for example, `x` and `X` are the same variable, so you need to pay attention here too. Also, the source has the tendency to be a little on the verbose side, but none of these are problems that should keep you from using it.

With all that said, FreeBASIC is an easy language, and the generated code is suprisingly fast.


## The Code
There are some modules in this folder. As with all programming languages, splitting the program in modules makes them easier to develop and mantain. I tried to follow the pattern in this implementation.

1. `complex.bas`: this is a module that implements a small complex library, since it is not a native data type in FreeBASIC. Modules in the language need to use `#ifndef` guards, inherited from C, to avoid multiple inclusions.

2. `test_it.bas`: a small module with functions to test and time the execution of the Fourier transforms. Correction checking is mostly done visually.

3. `fft.bas`: implements the Fourier Transforms. This module has a number of functions. `DirectFT` computes the Fourier transform using the definition, which is slow, but can be used for vectors of any length. `RecursiveFFT` and `IterativeFFT` implement the recursive and iterative Cooley-Tukey algorithm for vectors of power of 2 length (that is, 2, 4, 8, 16, 32, 64 and so on). `RecursiveNFFT` implements the Cooley-Tukey algorithm for vectors which length is a composite number.

4. `main_fft.bas`: main program to run power of 2 length FFTs a number of times and print the average running time for all implementations.

5. `main_anyfft.bas`: main program to run composite number length FFTs a number of times and print the average running time for all implementations.

## Compiling and Running
You will need to download FreeBASIC from their website if it is not installed in your machine. In my Ubuntu Linux box, it is not available from the package manager, so I needed to install it from external sources. The installation, however, is easy and you should have no problems with it.

Once the compiler is installed, you can compile the programs with the command:

```
$ fbc main_fft.bas
```

This will generate a `main_fft` executable in the same directory. Just run it and the program executes. To compile `main_anyfft.bas`, just substitute the name of the program in the same command line as above.
