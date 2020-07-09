# C++ Version

This folder contains the C++ implementation of the Discrete Fourier Transform. C++ was created in 1985 by Bjarne Stroustrup as an object-oriented extension to C. In its first years, it was very clumsy, mainly because C, the original language, is a very low level language, and object orientation is a high level concept. I remember once reading that C++ was trying to build an octopus by putting four tentacles in a dog. It was kind of true then, but I don't think it's like that anymore. The language, its libraries and community have grown, and the C++ environment is very mature and stable these days. It gained a lot of momentum when, back in the end of last century, it was seen as the perfect language to develop desktop applications.

Originally, C++ didn't have automatic memory management, and that is still true for primary data types and some classes derived from them. However, a lot of the classes in the standard library can manage their own memory, which makes it easier to deal with them. Still, this puts some distance from C++ to C, where everything must be managed by the programmer. Then again, it might not be the intent of C++ developers to follow the original philosophy. C++ stands on its own.

As happens with C, there are a lot of different implementations of the language. If you are working in a Linux system, your system most likely has `g++`, the GNU C++ Compiler, installed, and you can use it to compile and runs the programs on this folder; check your package manager. If you are in Windows or MacOS, however, you will probably have to install a development package. The companies developing these systems have them, choose whichever you like best.


## Comments on the Language

Developing with C++ poses a great difficulty: if you are used to program in C, it is very easy to *``fallback''* to C syntax and programming paradigm. If you don't know or forget how to do some simple task, C functions pop to your mind, and that might be a problem with your program. While this may sound silly, it *can* actually be a problem: there are a lot of things in both languages that look and sound similar, but are internally managed in completely different ways. This might lead to some incompatibilities that can make your program not work. Also, your program will end up as a frankenstein collation of low level C and object-oriented C++ trying to work toghether. It might not fail, but it will look ugly -- I tried to avoid this, but I can't be sure.

As a positive characteristic, C++ allows for operator overloading. This is a small characteristic that might not be that important or might be overused, but it does makes sense in programs that deal with mathematical objects, such as complex numbers or arrays.


## The Programs

There are two programs in this folder:

1. `fft.cpp`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.ppc`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. The Standard C++ library already have this, but I wanted to implement my own (that helps me to understand what the language can do). Also, if I was to follow the general guidelines, my complex library should come in a separate module, with a header file and so on. That would be extremely easy to do, but since these are very simple programs, I didn't think I needed that (I might change my mind in the future, however).


## Compiling and Running

Compile and run these programs will depend a lot on the compiler that you have installed. Since my system is Linux with the GNU C++ Compiler (`g++`) installed, I'll go with it. It is very easy to compile, just issue this command in the shell:

```
$ g++ -o fft fft.cpp -lm
```

to compile the `fft.cpp` file (don't forget the `-lm` switch to link the math library). This will generate a `fft` executable file in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run the `anyfft.cpp` file, follow the same steps, just change `fft` to `anyfft` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods.
