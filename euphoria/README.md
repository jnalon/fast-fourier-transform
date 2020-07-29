# Euphoria Version

This folder contains the Euphoria version of the Discrete Fourier Transform. Euphoria is a very simple imperative language created in 1993 for the Atari ST computer, and then convert to DOS. The main points of the design philosophy is that it should be simple, legible, performant and allow for rapid development. You can find more information about the language in its [website](https://openeuphoria.org/).

The language is indeed simple and legible -- in fact, it is designed in such a way that complex data structures, algorithms or libraries are responsability of the programmer. As an example, it has only two data types: `atom`, describing objects that can't be broken (such as numbers or characters) and `sequence`, lists of atoms (but you can also combine them to make your own data structure). This makes the core of the language very small and easy to learn and also very fast. I don't have experience developing larger software using Euphoria, but I suspect that this flexibility can be an issue if the code is not fully documented or designed with care beforehand.

Anyway, as you can see from the code below, Euphoria is a nice language -- I enjoyed writing code for it, though it took me some time to figure out one or other thing. Unfortunatelly, it didn't become very popular: a quick search on GitHub shows that the only project using Euphoria is the language development itself. But, if you want to learn a new language and already know the bigger ones, Euphoria might be a nice option.


## Comments on the Language

As said before, the main characteristic of Euphoria is its simplicity: two basic data types, so you have to do most of the work. This means that anything that is not trivial is responsability of the programmer. The language, however, provides constructs to allow you to build your own types in a very simple way. The type definition works different from other languages: while in traditional languages a data structure defines how the memory is used, in Euphoria it is a kind of a function that returns true or false, if the object follows or not the desired structured. That is: Euphoria defines *type checking* instead of *type structure*.

The resulting code is very readable -- in fact, I think only Python generates a more readable code (depending on what features of the language you use). Unfortunately, since most of the abstractions, such as structures, classes, operator overloading, and others, are not in the core of the language, this means that your code will probably have parts that define the guts of your objects that may look *ugly*. But the part of the code where you use them will be quite straightforward if you do your work well.

Last, the language is fast. It can be compiled or interpreted, but I got fast results even with the interpreted execution of the code. The compiler doesn't generate native code -- instead, it translates the Euphoria program to C, which can be compiled with any C compiler (such as `gcc`). I wasn't able to compile these programs, because some libraries are missing and I couldn't find in the documentation what development files should be installed in my system.



## The Programs

There are two programs in this folder:

1. `fft.ex`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.ex`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. I couldn't find a module for complex numbers in the standard library, but that wasn't a problem, since it is easy to define my own. It is unfortunate that Euphoria can't overload operators, however, since that would really improve readability in cases like this.


## Running

Euphoria is probably not installed in your system, and probably is not on your package manager too (in my Ubuntu system, however, I could install it directly from the package system). Installing it, however, is not complicated. Go to the home page of the language, download the program and follow the instructions and you're done.

To run the programs with the interpreter, you just issue the command:

```
$ eui fft.ex
```

The program will start and run. To run `anyfft.ex`, just change the `fft.ex` to `anyfft.ex` where needed.
