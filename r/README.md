# R Version

This folder contains the R version of the Discrete Fourier Transform. It was created in 1993 to deal with general numeric and statistical computations, and has found a great deal of success among staticians and data analysts. R puts emphasis on the development of numerical computations, science, plotting and other tools that help mathematical development. It is not, thus, a general purpose language, though there are modules to deal with things other than math.

It isn't a very powerful language, nor very efficient. Its power lies in the fact that the community is very active, and that new advanced functionality are promptly implemented and made available -- usually, by the scientists that developed the techniques themselves. This makes it very up-to-date and attractive to those dealing in data mining, data analysis, and machine learning in general.

To overcome the problems with efficiency, R can nativelly interface with C and Fortran, and with other languages as well through specific packages. Thus, a problem that runs slow in R can be converted to C and called without effort. Beside that, R doesn't have a lot of surprises or support for advanced programming techniques. It is a language created with an specific objective - model simulation -, and does that quite well.


## Comments on the Language

As any scripting language, R is very easy to learn. It has a simple syntax and a straight forward execution model. Among its features, R is capable of vectorizing operations, which allows to develop functions more by expressing *what* they do, instead of *how* they do. This is a characteristic of functional programming languages, and because of this, I would classify R as a *"semi-functional"* language.

On the other hand, R can be very slow. Although it is possible to enhance its performance by writing parts of code in C or Fortran (as interfacing with those languages is easy), R native code is not very efficient. Nevertheless, much of the code in R libraries is written in R itself. Overall, R is an almost perfect choice if you want to deal with heavy computation -- there are faster languages, however.

I didn't have any surprises writing the code, and people used to any imperative programming language will be at home with R, only having to learn its syntax. It is well possible that the language has some unique characteristics that I didn't use that can make it more efficient or easy to write.


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
