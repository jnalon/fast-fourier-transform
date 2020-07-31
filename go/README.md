# Go Version

This folder contains the Go version of the Discrete Fourier Transform. Go was created by Google in 2009. It was created to improve development at Google, and retains most of the characteristics that are usually found in compiled programming languages: it is staticly typed, it runs efficiently, but is also very readable, and supports multiprocessing. The most striking feature is that, contrary to the most modern languages developed in the last 10 years, Go doesn't have classes -- although object orientation is supported in other ways.

The language is also very easy to learn, at least to start producing functioning code. You probably need some experience to deal with the most complicated aspects of it -- but every programming language works this way. An interesting aspect of the language is the support for multiprocessing, but I can't comment a lot on it, since I didn't use it. You can find more information at their [website](https://golang.org/).


## Comments on the Language

Go is easy to learn and to write code. Its syntax is not too different from what a seasoned programmer is already used to, but they diverge at some point. I found array declaration kind of weird, with squared brackets *before* type name, which is unusual and, frankly, looked like Go designers were only trying to be different. The assignment/declaration operator (`:=`) is welcome, since you don't need to clog your code with type declarations. Aside from that, there isn't a lot of surprises in the way code is written. I found, however, that it lack some abstractions that would make it wonderful, such as initializers for vectors or a better syntax for slices. Also, Go doesn't promote *up* numeric types (*ie*, convert integer to float when needed), so you have to add a lot of type casting. That makes the code look ugly at some points.

There are full support for creation of modules and libraries -- I, in fact, wrote a small library for image processing, but it wasn't that easy to find out the better way to structure the source files in directory structure. The generated code is very fast and, although I didn't test it, is probably memory efficient also. I imagine that there is full support for in-code documentation.


## The Programs

There are two programs in this folder:

1. `fft.go`: this implements `DirectFT`, `RecursiveFFT` and `IterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.go`: this implements `DirectFT` and `RecursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `DirectFT`, and shows no gain in efficiency at all.

Go has a complex number library. This is fortunate because, since the language doesn't support operator overloading, complex arithmetic code would probably look very cumbersome.


## Running

If you use Linux, Go might already be installed on your system, or can be easily installed via package manager. In Windows or MacOS, you will have to download, install and configure the packages.

Go is compiled, but you can compile and run the programs with only one command:

```
$ go run fft.go
```

The program will start and run. To run `anyfft.go`, just change `fft.go` to `anyfft.go` where needed.
