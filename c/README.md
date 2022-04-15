# C Version

This folder contains the implementation of the Fast Fourier Transform using the C Programming Language. C is a general purpose procedural structured programming language. It was created by Dennis Ritchie in 1972 with the explicit purpose of developing the Unix Operating System at Bell Labs. It has found a great deal of success since then, and most of present day languages, if not all of them, are heavy influenced by it.

C was designed to be a structured language working *very* close to the machine, and there are people who describe it as a *high-level assembly*. Programming in the language allows you to direct access memory contents and devices with a very simple syntax, which makes it very flexible; C compilers are often heavily optimized, so the code is small and runs fast. Those things, however, come with a price: modern capabilities such as automatic memory management or data structures are not present, they have to be programmed or included with external libraries. This makes it very difficult to debug C code, and the amount of side effects -- especially when dealing with pointers -- can make C programs very unstable. However, there is a huge number of external and third party libraries that help deal with the problem (which can lead to other problems, such as incompatibility between libraries, and so on, but that's a discussion beyond the scope of this document).


## Comments on the Language

First, I must say that I've been using C most of my life, so it was exactly difficult for me. I had to check here and there for help with a function of library, but that was all. However, I don't claim that this is the most efficient implementation (it isn't, of course), but I explicitly tried to avoid some tricks that experienced programmers use. I think that resulted in a very readable and easy to understand code.

Strictly speaking, C is *not* an easy language. Although you can program in it using the basic structures, you won't harvest the full power of the language - which is probably what happened with my code here. There are a lot of techniques and idioms that will make the programs *safer* and *faster*. Knowing and learning such techniques demands a lot time spent on studying and practicing with the language, and not all programmers will agree on what the best tricks are, or whether they are good to use or not.

However, since most languages are at least in part based on it, C won't feel different, even if it was the first time using it (although, if that is the case, I'd recommend starting with easier algorithms). In fact, C won't feel that complicated if you avoid problematic things like memory management. Every program, however, *will* have to deal with memory with some extent. There are two techniques that I use to mitigate the problems: first, every variable that needs it has its memory allocated at the beginning of the scope and is freed in the end of the scope, so it's easy to figure if something is not or wrongly allocated or deallocated; second, if a function must return an array or a set of values, the memory is allocated and freed outside of the function. These are, in fact, good practices in C.


## The Programs

There are four programs in this folder:

1. `fft.c`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.c`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

3. `fft_complex.c`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.). The difference of this program from the first one is that this uses the native complex library of the C language.

2. `anyfft_complex.c`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all. The difference of this program from the first one is that this uses the native complex library of the C language.

Besides the transform functions, the first two files also implement a small library to deal with complex numbers. The Standard C library already have a `complex.h` module, but I wanted to implement my own (that helps me to understand what the language can do). Also, if I was to follow the general guidelines, my complex library should come in a separate module, with a header file and so on. That would be extremely easy to do, but since these are very simple programs, I didn't think I needed that (I might change my mind in the future, however). The other two use the Standard C `complex.h` library, but the differences in performance were negligible, if there's any at all.

The good practices of development in C also mandate that I should actually build a header file (`.h`) to hold the functions and include that in the main program. The header file should control the size of data and precision of operations. That is correct, of course, but if I were to do that, I would be diverging from my main intent, that was creating the Fast Fourier Transform code. The functions, anyway, can be transfered to bigger projects with a more adequate structure.


## Compiling and Running

Compiling and running these programs will depend a lot on the compiler that you have installed. Since C is an open language, there are numerous different compilers. If you are using a Linux OS, chances are that you have already a C compiler installed, most likely `gcc`; check in the software manager of your system. If you are using Windows or MacOS, however, you will need to install a development package. In these cases, you can choose between the open ones (such as `gcc` from [MinGW](http://mingw.org/) or the development packages made available by the producer themselves. Choose the one you like most.

My system is Linux with the GNU C Compiler (`gcc`) installed, I'll go with it. It is very easy to compile, just issue this command in the shell:

```
$ gcc -o fft fft.c -lm
```

to compile the `fft.c` file (don't forget the `-lm` switch to link the math library). This will generate an executable file named `fft` in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run `anyfft.c`, `fft_complex.c` or `anyfft_complex.c`, follow the same steps, just change the program name in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods.
