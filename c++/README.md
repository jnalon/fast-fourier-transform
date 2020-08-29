# C++ Version

This folder contains the C++ implementation of the Discrete Fourier Transform. It is a general purpose procedural and object oriented language, created in 1985 by Bjarne Stroustrup as an object-oriented extension to C, and because of that, most of C programs *could* be compiled by C++ compilers. However, since C is a low-level language, and object oriented is a high-leve concept, it used to be very clumsy, and programs developed in the language were riddled with memory leaks and other hard to find bugs. I remember once reading that C++ was trying to build an octopus by putting four tentacles in a dog.

It was true then, but I don't think it's like that anymore. It gained a lot of momentum when, back in the end of last century, it was seen as the perfect language to develop desktop applications. The language, its libraries and community have grown, and the C++ environment is very mature and stable these days.

Originally, C++ didn't have automatic memory management, and that is still true for primary data types and some data structures and classes. However, a lot of the classes in the standard library can manage their own memory, which makes it easier to deal with them. This puts some distance from C++ to C, where everything must be managed by the programmer, but, probably, it isn't a concern of C++ developers to follow the original philosophy anymore. C++ stands on its own.


## Comments on the Language

Developing with C++ poses a difficulty (that happens with any derived language): if you are used to program in C, it is very easy to *"fallback"* to C programming, and you will be *"programming in C, using C++ syntax"*. If you don't know or forget how to do some simple task, you will feel tempted to use C functions and usual solutions. While this may sound silly, it *can* actually be a problem: there are a lot of things in both languages that look and sound similar, but are internally managed in completely different ways. This can lead to incosistencies that make your program not work. Also, your program will end up as a *"frankenstein"* collation of low level C and object-oriented C++ trying to work toghether. It will most likely not fail, but it will look ugly.

C++ allows operator overloading, multiple inheritance and other features that might be considered *dangerous*. While that is true, a good design can avoid most if not all of the problems. In the development of these programs, multiple inheritance wasn't necessary, but operator overloading was very helpful. Without this, an expressions such as

```
X[k] = X[k] + Wkn * X[k]
```

where all variables are complex, would have to be written as:

```
X[k] = X[k].add(Wkn.mul(X[k]))
```

which can get ugly for more complex expressions.


## The Programs

There are two programs in this folder:

1. `fft.cpp`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.ppc`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. The Standard C++ library already have this, but I wanted to implement my own (that helps me to understand what the language can do). Also, if I was to follow the general guidelines, my complex library should come in a separate module, with a header file and so on. That would be extremely easy to do, but since these are very simple programs, I didn't think I needed that (I might change my mind in the future, however).

As a last note, the programs have *few* characteristics of object orientation. This is because the Fast Fourier Transform is better implemented as an operation (and, thus, as a function) than as a method of a class. In fact, to do it in that way, I would have to create a class to hold the vector data and implement some additional methods to create, allocate and dispose memory and so on. While I could have done this, that would diverge from my first intent, that was to implement the Fast Fourier Transform. So, you might argue that this is - as I said above - C written with C++ syntax, but the functions can be easily transfered to bigger class oriented projects.


## Compiling and Running

As happens with C, there are a lot of different implementations of the language. If you are working in a Linux system, your system most likely has `g++`, the GNU C++ Compiler, installed, and you can use it to compile and run the programs on this folder; just check the package manager for your distribution, and install the needed packages. If you are in Windows or MacOS, however, you will probably have to install a development package. There are a lot of options (such as [MinGW](http://mingw.org/) for Windows, and GNU C++ is probably available for MacOS. As an alternative, Windows has [Visual Studio](https://visualstudio.microsoft.com/) and MacOS has [XCode](https://developer.apple.com/), both of them with free versions that you can download and use.

I'm in a Linux system, so I will use G++. To compile the program, just issue the command:

```
$ g++ -o fft fft.cpp -lm
```

to compile the `fft.cpp` file (don't forget the `-lm` switch to link the math library). This will generate an executable file named `fft` in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run the `anyfft.cpp` file, follow the same steps, just change `fft` to `anyfft` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods.
