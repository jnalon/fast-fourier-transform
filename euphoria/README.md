# Euphoria Version

This folder contains the Euphoria version of the Discrete Fourier Transform. Euphoria is a very simple procedural language created by Robert Craing in 1993 for the Atari ST computer, and then converted to DOS. The main design goal was that it should be simple, legible, performant and allow for rapid development. In fact, one of its main points is that *it shouldn't provide* any complicated abstraction that might reduce its efficiency. You can find more information about the language in its [website](https://openeuphoria.org/).

Euphoria is very flexible, as it can be compiled to bytecode, translated to C, or interpreted. Executables can be made either by compilation of by binding the interpreter and runtime toghether, which you can use to make Euphoria a nice scripting language for bigger programs. Unfortunatelly, it didn't become very popular: a quick search on [GitHub](https://github.com/) shows that the only project using Euphoria is the language development itself. But, if you want to learn a new language and already know the bigger ones, Euphoria might be a nice option.


## Comments on the Language

The language is indeed simple and legible, it almost look like an script language. It has only two data types: `atom`, describing objects that can't be broken (such as numbers or characters) and `sequence`, lists of atoms (but you can also combine them to make your own data structure). This makes the core of the language very small and easy to learn and also very fast. This, however, has a drawback: whatever is not trivial is always responsability of the programmer, and even common data structures (such as queues or stacks) must be programmed by hand. I feel, however, that it *could* have an abstraction or two, such as namespaces or better module management. It might not be a problem for easy structures, but I suspect that this flexibility can be an issue if the code is not fully documented or designed with care beforehand.

The language provides constructs to allow you to build your own types in a very simple way. There is a difference from how type definition works on common languages from how it works here: while in traditional languages a data structure defines how the memory is used, in Euphoria it is a kind of a function that returns true or false, if the object follows or not the desired structured. That is: in Euphoria, you define *type checking* instead of *type structure*, and it is not mandatory - you can build functions that don't type check at all. I think this is a nice solution to the problem of allowing dynamic typing but must do type checking at runtime.

The resulting code is very readable -- in fact, I think only Python generates a more readable code (depending on what features of the language you use). Unfortunately, since most of the abstractions, such as structures, classes, operator overloading, namespace and others, are not in the core of the language, this means that your code will probably have parts that define the guts of your objects that may look *ugly*. But the part of the code where you use them will be quite straightforward if you do your work well.

You can *emulate* most of those things, but it is *your* responsability. You can create something that works like namespaces by structuring your code in different files and following a naming convention; something similar could be done to emulate object orientation. For example, if you're dealing with image processing, you could create a module name (say) ``image.ex``, and prefix every function with ``img_``. The language, however, doesn't provide with a semantic to do that, nor help you manage it, so this can be error prone, and more difficult in the end.

Last, the language is fast. It can be compiled or interpreted, but I got fast results even with the interpreted execution of the code. The compiler doesn't generate native code -- instead, it translates the Euphoria program to C, which can be compiled with any C compiler (such as `gcc`). I wasn't able to compile these programs, because some libraries are missing and I couldn't find in the documentation what development files should be installed in my system.


## The Programs

There are two programs in this folder:

1. `fft.ex`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.ex`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. I couldn't find a module for complex numbers in the standard library, but that wasn't a problem, since it was easy to define my own. It is unfortunate that Euphoria can't overload operators, however, since that would really improve readability in cases like this.


## Running

If you use Linux, Euphoria is probably not installed in your system, and there is a good chance that it isn't listed on your package manager. It is listed as one of available packages in Ubuntu, however. Windows and MacOS users will have to download the packages and perform the installation, but, in any case, this is not a problem  , and probably is not on your package manager too (in my Ubuntu system, however, I could install it directly from the package system). Installing it, however, is not complicated. Go to the home page of the language, download the program and follow the instructions and you're done.

To run the programs with the interpreter, you just issue the command:

```
$ eui fft.ex
```

The program will start and run. To run `anyfft.ex`, just change the `fft.ex` to `anyfft.ex` where needed. I wasn't able to translate to C - the translation process finished without any problems, but I was lacking some header files that I couldn't find. This is *most likely* because Euphoria is trying to compile for 32-bits in a 64-bits operating system. I will probably look further in the problem soon.
