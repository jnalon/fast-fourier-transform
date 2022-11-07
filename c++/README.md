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

One of the great resources in C++ is the use of templates. If you check the history of these implementations, you will see that the FT was implemented for `float`s only -- but it was extremely easy to adapt classes and functions to use templates. Thus, with minimal effort, I was able to implement the FT both for `float`s and `double`s. With the possibility of using operator overloading, I was tempted to try a fixed point implementation. I succumbed to the temptation, of course, and the result can be seen in the source code.

## The Programs
There is a number of files in this folder. The modularization, while makes it easier to develop, understand the program, can cause the adverse effect of multiplying the number of project files. But they are easy to follow when documented:

1. `my_complex.h`: this implements a very simple and small complex number library. This was implemented as a templated class, so you can choose the type of real and imaginary parts. This made going from `float` to `double` to `FixedPoint` incredibly easy. However, since the class is templated, declaration and definition must be in the same file, so there is no `.cpp` for this header.

2. `test_it.h`: implementation of functions used to measure the time spent in the computation of the FFT. Since these functions are also templated, declaration and definition can be found in the same file.

3. `fixed_point.h` and `fixed_point.cpp`: respectively, header and implementation of a small fixed point arithmetic library. This implements fixed point numbers as a class encapsulating a `long int` to hold the number values, and the main operations with these numbers. There are also functions to evaluate Taylor series, which are used to compute _sine_ and _cosine_ of fixed point numbers.

4. `fft.h`: this implements `direct_ft`, `recursive_fft`, `iterative_fft` and `recursive_nfft`. The first function computes the FFT for vectors of any length, the following two use Cooley-Tukey decomposition to compute for vectors of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.), and the last one use Cooley-Tukey decomposition to compute for vectors with length a composite number. As in the case of the first two files, the functions are templated, so declaration and definition can be found here.

5. `fp_test.c`: a small program to test the fixed point library by printing the results of some operations. The checking is mainly done visually.

6. `main_fft.c`: main program to run power of 2 length FFTs a number of times and print the average running time for all implementations of the FFT.

7. `main_anyfft.c`: main program to run composite number length FFTs a number of times and print the average running time for all implementations of the FFT.

8. `main_test.c`: this is a simple program to compute the FFT for vectors of some lengths, using all implementations available on the libraries, and printing the results to check if the results are right. A more correct way to implement this was to actually build a test suite, but that could get too complicated, and I don't think it's needed for a small and simple piece of code as this.

If you checked these programs earlier, you probably noticed that they weren't written in the best C++. While this was intentional to save memory and execution time, I felt that I needed to refactor the programs to a better use of the capabilities of the language. I think they're resembling better C++ code now.

## Compiling and Running
As happens with C, there are a lot of different implementations of the language. If you are working in a Linux system, your system most likely has `g++`, the GNU C++ Compiler, installed, and you can use it to compile and run the programs on this folder; just check the package manager for your distribution, and install the needed packages. If you are in Windows or MacOS, however, you will probably have to install a development package. There are a lot of options (such as [MinGW](http://mingw.org/) for Windows, and GNU C++ is probably available for MacOS). As an alternative, Windows has [Visual Studio](https://visualstudio.microsoft.com/) and MacOS has [XCode](https://developer.apple.com/), both of them with free versions that you can download and use.

I'm in a Linux system, so I will use G++. To compile the program, just issue the command:

```
$ g++ -o fft main_fft.cpp fixed_point.cpp -lm
```

to compile the `main_fft.cpp` file (don't forget the `-lm` switch to link the math library). This will generate an executable file named `fft` in the same folder, that can be run with the command:

```
$ ./fft
```

To compile and run the `main_anyfft.cpp` file, follow the same steps, just change `main_fft` to `main_anyfft` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods. The test program can be compiled by substituting `main_test.cpp` in the command line.

If you want to use this implementation in your code, it's easy: just copy the needed files. You might want to tweak the `my_complex.h` and `fft.h` to use a specific data type (I suggest using `double`, as most processors have instructions to deal with double precision floats). If you want to use the fixed point library, the procedure is the same.

**However**, with all that being said, I don't recomend using this code, as its purpose is only educational. You should use more stable and tested libraries, such as [FFTW](http://www.fftw.org/).

