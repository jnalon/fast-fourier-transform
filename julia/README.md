# Julia Version

This folder contains the Julia version of the Discrete Fourier Transform. Julia is dynamic programming language created in 2012, especializing in scientific computation. It is clearly inspired in Matlab and Python, but there are a lot of constructs that differ from both languages. The main feature, at least as I see, is a just-in-time compiler that aims to provide fast execution times. Also, it can nativelly call C and Fortran functions, and interface with a lot of other languages, such as Python and Java, through libraries. It is being used by a lot of business and research groups for their numerical computation, but its following still has to grow.

The language itself is very simple, and any person who already programmed in Matlab, Python or similar scripting languages will feel at home with Julia. Especially because the language provides some resources such as vectorization of operations, which leads to very compact, readable and efficient code. Julia, however, is not object-oriented in the strict sense, thoug there are abstractions that allow to deal with the paradigm. You can find more information on its [website](https://julialang.org/)


## Comments on the Language

Julia is easy to learn and write. However, because it is very similar to other scripting languages, programmers might find that they are, instead of using the better of Julia's capacities, they are programming in their older languages using Julia's syntax. That is, of course, a problem with every programming language, and I'm sure that there are better ways to write the code in these files.

The code, maybe because of this, didn't run as fast as I expected. It is very easy to learn a new language, but it might be very difficult use it in the best. But I found that even the direct FT, which implementation is straight forward, ran way slower than i hoped. One other problem that I found was that, at some expressions, vectorization didn't work quite well, or not at all. Again, it might be that I wasn't able to use the language's resources, but I was kind of disappointed.


## The Programs

There are two programs in this folder:

1. `fft.jl`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.jl`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

There was no need to write functions for complex numbers, as they're native in the language.


## Running

If you use Linux, Julia can be installed from your package manager; you might have to install a package or two from the language command line, but that is easy. On Windows or MacOS, you will have to download the packages, but installation is easy.

Julia is not interpreted, it is converted to intermediate code that is compiled to run. In any case, you can run the programs by just issuing the command:

```
$ julia fft.jl
```

To run `anyfft.jl`, just change `fft.jl` where needed.
