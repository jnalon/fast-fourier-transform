# C Version

This folder contains the implementation of the Fast Fourier Transform using the C Programming Language. C was created by Dennis Ritchie in 1972 with the explicit purpose of developing the Unix Operating System at Bell Labs. It has found a great deal of success since then, and most of present day languages, if not all of them, are heavy influenced by it.

C was designed to be a structured language working *very* close to the machine. there are people who describe it as a *high-level assembly*. Programming in the language allows you to direct access memory contents and devices with a very simple syntax, which makes it very flexible; C compilers are often heavily optimized, so the code is small and runs fast. Those things, however, come with a price: modern capabilities such as automatic memory management or data structures are not present, they have to be programmed. This makes it very difficult to debug C code, and the amount of side effects -- especially when dealing with pointers -- can make C programs very unstable. However, there is a huge number of external and third party libraries that help deal with the problem (which can lead to other problems, such as incompatibility between libraries, and so on, but that's a discussion beyond the scope of this document).

Since C is an open language, there are numerous different compilers. If you are using a Linux OS, chances are that you have already a C compiler installed, most likely `gcc`; check in the software manager of your system. If you are using Windows or MacOS, however, you will need to install a development package. Choose the one you like most.


## The Programs

There are two programs in this folder:

1. `fft.c`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.c`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. The Standard C library already have a `complex.h` library, but I wanted to implement my own (that helps me to understand what the language can do). Also, if I was to follow the general guidelines, my complex library should come in a separate module, with a header file and so on. That would be extremely easy to do, but since these are very simple programs, I didn't think I needed that (I might change my mind in the future, however).


## Compiling and Running

Compile and run these programs will depend a lot on the compiler that you have installed. Since my system is Linux with the GNU C Compiler (`gcc`) installed, I'll go with it. It is very easy to compile, just issue this command in the shell:

```
$ gcc -o fft fft.c -lm
```

to compile the `fft.c` file (don't forget the `-lm` switch to link the math library). This will generate a `fft` executable file in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run the `anyfft.c` file, follow the same steps, just change `fft` to `anyfft` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods.


## Comments on the Language

First, I must say that I've been using C most of my life, so it was exactly difficult for me. I had to check here and there for help with a function of library, but that was all. However, I don't claim that this is the most efficient implementation (it isn't, of course), but I explicitly tried to avoid some tricks that experienced programmers use. I think that resulted in a very readable and easy to understand code.

Strictly speaking, C is *not* an easy language. However, since most languages are at least in part based on it, it won't feel different, even if it was the first time using it (although, if that is the case, I'd recommend starting with easier algorithms). In fact, C won't feel that complicated if you avoid problematic things like memory management. Every program, however, *will* have to deal with memory with some extent. There are two techniques that I use to mitigate the problems: first, every variable that needs it has its memory allocated at the beginning of the scope and is freed in the end of the scope, so it's easy to figure if something is not or wrongly allocated or deallocated; second, if a function must return an array or a set of values, the memory is allocated and freed outside of the function. These are, in fact, good practices in C.

