# C Version
This folder contains the implementation of the Fast Fourier Transform using the C Programming Language. C is a general purpose procedural structured programming language. It was created by Dennis Ritchie in 1972 with the explicit purpose of developing the Unix Operating System at Bell Labs. It has found a great deal of success since then, and most of present day languages, if not all of them, are heavy influenced by it.

C was designed to be a structured language working *very* close to the machine, and there are people who describe it as a *high-level assembly*. Programming in the language allows you to direct access memory contents and devices with a very simple syntax, which makes it very flexible; C compilers are often heavily optimized, so the code is small and runs fast. Those things, however, come with a price: modern capabilities such as automatic memory management or data structures are not present, they have to be programmed or included with external libraries. This makes it very difficult to debug C code, and the amount of side effects -- especially when dealing with pointers -- can make C programs very unstable. However, there is a huge number of external and third party libraries that help deal with the problem (which can lead to other problems, such as incompatibility between libraries, and so on, but that's a discussion beyond the scope of this document).

## Comments on the Language
First, I must say that I've been using C most of my life, so it wasn't exactly difficult for me. I had to check here and there for help with a function of library, but that was all. However, I don't claim that this is the most efficient implementation (it isn't, of course), but I explicitly tried to avoid some tricks that experienced programmers use. I think that resulted in a very readable and easy to understand code.

Strictly speaking, C is *not* an easy language. Although you can program in it using the basic structures, you won't harvest the full power of the language - which is probably what happened with my code here. There are a lot of techniques and idioms that will make the programs *safer* and *faster*. Knowing and learning such techniques demands a lot time spent on studying and practicing with the language, and not all programmers will agree on what the best tricks are, or whether they are good to use or not.

However, since most languages are at least in part based on it, C won't feel different, even if it was the first time using it (although, if that is the case, I'd recommend starting with easier algorithms). In fact, C won't feel that complicated if you avoid problematic things like memory management. Every program, however, *will* have to deal with memory with some extent. There are two techniques that I use to mitigate the problems: first, every variable that needs it has its memory allocated at the beginning of the scope and is freed in the end of the scope, so it's easy to figure if something is not or wrongly allocated or deallocated; second, if a function must return an array or a set of values, the memory is allocated and freed outside of the function. These are, in fact, good practices in C.

Previous versions of the files in these directories (check the history) had separated implementations using a small custom complex number library and using C native complex numbers. But, since it is a good idea to compare the results of both implementations, it was refactored to a more modularized version that allowed to mix the functions. This lead to a problem: since C doesn't have _namespaces_ (or something equivalent to it), function names needed to be prefixed to be different. This lead to some code duplication, which is _never_ a good thing. Here, it isn't a big problem: if you are using one of the implementations of the FFT, you can disregard the other, and adapt the names as you wish. _There are_ ways to solve this, but that would lead to very complex code, and that is not needed in a project so small as this.

## The Code
There are some files in this folder. A good practice in C language is to modularize your program by separating uncorrelated functions. While this may seem that a lot of unnecessary files are created, it is easier to mantain and find bugs, and also to copy only the needed files to a project. For example, if you only need the FFT using native C complex numbers, only the corresponding file needs to be added to your project:

1. `my_complex.h` and `my_complex.c`: respectively, the header and implementation files for a very simple and small complex number library used to compute the FFT.

2. `test_it.h` and `test_it.c`: respectively, the header and implementation files for functions used to measure the time spent in the computation of the FFT. In this file, there are test functions using the `my_complex` library defined in the files described above, and functions using C native complex numbers with single and double precision.

3. `fft.h` and `fft.c`: respectively, the header and implementation files for the functions computing the FFT using the small `my_complex` library. These files implements `direct_ft`, `recursive_fft`, `iterative_fft` and `recursive_nfft`. The first function computes the FFT for vectors of any length, the following two use Cooley-Tukey decomposition to compute for vectors of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.), and the last one use Cooley-Tukey decomposition to compute for vectors with length a composite number.

4. `fft_float.h` and `fft_float.c`: respectively, the header and implementation files for the functions computing the FFT using C language native complex numbers. These files implements `float_direct_ft`, `float_recursive_fft`, `float_iterative_fft` and `float_recursive_nfft`. The first function computes the FFT for vectors of any length, the following two use Cooley-Tukey decomposition to compute for vectors of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.), and the last one use Cooley-Tukey decomposition to compute for vectors with length a composite number.

5. `fft_double.h` and `fft_double.c`: respectively, the header and implementation files for the functions computing the FFT using C language native complex numbers with double precision. These files implements `double_direct_ft`, `double_recursive_fft`, `double_iterative_fft` and `double_recursive_nfft`. The first function computes the FFT for vectors of any length, the following two use Cooley-Tukey decomposition to compute for vectors of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.), and the last one use Cooley-Tukey decomposition to compute for vectors with length a composite number.

6. `main_fft.c`: main program to run power of 2 length FFTs a number of times and print the average running time for all implementations of the FFT.

7. `main_anyfft.c`: main program to run composite number length FFTs a number of times and print the average running time for all implementations of the FFT.

8. `main_test.c`: this is a simple program to compute the FFT for vectors of some lengths, using all implementations available on the libraries, and printing the results to check if the results are right. A more correct way to implement this was to actually build a test suite, but that could get too complicated, and I don't think it's needed for a small and simple piece of code as this.

## Compiling and Running
Compiling and running these programs will depend a lot on the compiler that you have installed. Since C is an open language, there are numerous different compilers. If you are using a Linux OS, chances are that you have already a C compiler installed, most likely `gcc`; check in the software manager of your system. If you are using Windows or MacOS, however, you will need to install a development package. In these cases, you can choose between the open ones (such as `gcc` from [MinGW](http://mingw.org/)) or the development packages made available by the manufacturer themselves. Choose the one you like most.

My system is Linux with the GNU C Compiler (`gcc`) installed, I'll go with it. It is very easy to compile, just issue this command in the shell:

```
$ gcc -o fft main_fft.c fft.c fft_float.c fft_double.c my_complex.c test_it.c -lm
```

to compile the `main_fft.c` file (don't forget the `-lm` switch to link the math library). This will generate an executable file named `fft` in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run `main_anyfft.c`, follow the same steps, just change `main_fft.c` to `main_anyfft.c` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods. You won't be surprised to find out that the native implementation of complex numbers is faster (but by a _very small_ margin). The test program can be built by substituting `main_test.c` for `main_fft.c` in the command line.

If you want to use this implementation in your code, it's easy. You probably will like to use the implementation on `fft_float.*` or `fft_double.*`. I **strongly** suggest that you test these libraries, to determine which of them runs faster on your system. I was happily surprised to see that the _double_ version was more performant -- that is probably because the processors in my machine have instructions to deal with double precision floats. Since my computer is very standard, this is probably your case too. Check the comments on each file and adapt your program to use them. You can use the test functions in `test_it.c` as a template on how to use them.

**However**, with all that being said, I don't recomend using this code, as its purpose is only educational. You should use more stable and tested libraries, such as [FFTW](http://www.fftw.org/).
