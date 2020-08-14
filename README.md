# Fast Fourier Transform

## Simple implementation of the Fast Fourier Transform in various languages

This project was created as a personal register of my effort to learn programming languages. When I'm learning a new language, I try a certain number of algorithms to see how much I was able to understand, and how much time it takes until I can produce code. Among the algorithms I always try are: Fibonacci sequence, Factorial, Quicksort, The Knight's Round, Towers of Hanoi, Fast Fourier Transform, and others. Since I'm an Electronic Engineer with heavy training in Computers and Telecomunications, I chose the FFT to investigate with more details.

Every implementation here does basically the same: they execute the Cooley-Tukey decomposition FFT for vectors of various lengths, almost always a power of two, although some implementations can deal with composite numbers. The transforms are run a number of times, and execution time is measured and averaged. A table with a comparison of efficiency is printed on screen.

Currently, this repository contains implementations in the languages on the list below. There is a directory for each language, so you can clone only what you want to see. Some folders have more than one version of a file -- usually, an implementation of the Cooley-Tukey algorithm for composite numbers (but not all directories have them).

Instructions on how to compile and run are in the comments of each file. There is no `make` file because those are very simple commands. Some of the implementations have dependencies, although I tried to keep them at a minimum. Again, those are in the comments in the header of each file.

I don't claim that any of these implementations is the most efficient way to deal with the FFT in any language. In fact, most likely, it isn't, as most of them have bindings for a tested and fast library, such as FFTW.


### Imperative Languages

The languages below are imperative languages. It's the most common programming paradigm, and most of the languages in existence are designed and implemented to support it. They're caracterized by sequences of commands, mutable variables, repetition and decision comands, and strict evaluation. Most of these languages have constructs to support a particular aspect of programming, but the implementation in most of them are very simmilar. I included here object oriented languages, since object orientation is just one other way to structure imperative programs:

* C
* C#
* C++
* Euphoria
* Go
* Groovy
* Java
* Kotlin
* Pascal
* Ruby


### Semi-functional Languages

I'm calling *semi-functional* those languages that, although they follow the imperative paradigm, have constructs that allow the programmer to think not in a sequence of commands, but in the relations between values. They still have the same characteristics of the languages above, but have some constructs that allow you to think in terms of functions instead of procedures:

* Julia
* MATLAB/Octave
* Python
* R
* Scala


### Functional Languages

Functional languages support the paradigm that programs should be thought of and designed as mathematical equations. Usually, variables are not mutable, there are no loop and decision commands. The supporters of the paradigm *claim* that programs designed with such philosophy are not prone to errors -- but, of course, that's not true. They also claim that functional programming is *"purer"* -- which is obvious also not true. Programming in a functional paradigm, however, is very different from the imperative paradigm, and it takes some time to get used, but it is something that every developer should know:

* Erlang


----
José Alexandre Nalon
