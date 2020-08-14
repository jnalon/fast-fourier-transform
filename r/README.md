# R Version

This folder contains the R version of the Discrete Fourier Transform. It was created in 1993 to deal with statistical computations. It is being widely used in data mining and analysis. It puts emphasis on the development of numerical computations, science, plotting and other tools that help mathematical development. It is not, thus, a general purpose language, though there are modules to deal with things other than math.

It isn't a very powerful language, nor very efficient. Its power lies in the fact that the community is very active, and that modules with new development in statistics are promptly implemented and made available -- normally, by the scientists that developed the techniques themselves. This makes it very attractive to those dealing in the field.


## Comments on the Language

As any scripting language, R is very easy to learn. It has a simple syntax and a straight forward execution model. Among its features, R is capable of vectorizing operations, which allows to develop functions more by expressing *what* they do, instead of *how* they do. This is a characteristic of functional programming languages, and because of this, I would classify R as a *"semi-functional"* language.

On the other hand, R can be very slow. Although it is possible to enhance its performance by writing parts of code in C or Fortran (as interfacing with those languages is easy), R native code is not very efficient. Nevertheless, much of the code in R libraries is written in R itself. Overall, R is an almost perfect choice if you want to deal with heavy computation -- there are faster languages, however.


## The Programs

There are two programs in the folder:

1. `fft.r`: this program implements `direct_ft`, `recursive_fft` and `iterative_fft` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms, and also with the internal implementation of the FFT available in the language. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.r`: this program implements `direct_ft` and `recursive_fft`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `direct_ft`, and shows no gain in efficiency at all.

This program didn't implement a complex arithmetic library, because complex numbers are native to the language.


## Running

R is an interpreted language. To run the program, type the command:

```
$ Rscript fft.r
```

To run `anyfft.r`, just substitute the filename on the command line.
