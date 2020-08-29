# Julia Version

This folder contains the Julia version of the Discrete Fourier Transform. Julia is dynamic programming language created in 2012, specializing in scientific computation. It is clearly inspired in Matlab and Python, but there are a lot of constructs that differ from both languages. It can nativelly call C and Fortran functions, and interface with a lot of other languages, such as Python and Java, through libraries. It is being used by a lot of business and research groups for their numerical computation, but its following still has to grow.

The language's [website](https://julialang.org/) says that Julia is actually multipurpose. That might be true, but it is very easy to see that it is slanted in the direction of scientific computing. Among its features, Julia has dynamic dispatching (although not in an object orientation context -- this is not supported by Julia), a built-in package manager and nice metaprogramming features. Also deserving mention is its support for coroutines. For ease of development, Julia has an in-built package manager and documentation system (much like Python's docstrings).

It is very simple, and any person who already programmed in Matlab, Python or similar scripting languages will feel at home with Julia. But it can be very powerful; for example, one of its features (that was usefull in these programs) was vectorization of operations, which leads to very compact, readable and efficient code. Julia also sports some support for functional programmng. You can find more information on its [website](https://julialang.org/)


## Comments on the Language

Julia is easy to learn and write. However, because it is very similar to other scripting languages, programmers might find that they are, instead of using the better of Julia's capacities, they are programming in their older languages using Julia's syntax. I'm sure, however, that the intention of the designers was not be different for the sake of being different: it is clear from using the language that what is proven to work (such as an *if* command), is written and works exactly as expected.

The execution of the programs is fast, but not as fast as you would hope from a language that has a just-in-time compiler. It is enough to get things done in an efficient way, but I was hoping it was faster. However, there probably are ways to optimize the code and make it run in a shorter length of time.


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
