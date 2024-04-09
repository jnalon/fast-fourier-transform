# Java Version
This folder contains the Java version of the Discrete Fourier Transform. Java is an object oriented language created by James Gosling at Sun Microsystems in 1995. It was designed with the goal of allow developing portable applications -- it's *motto* was *\"write once, run everywhere\"*. Although it inherited a lot of syntax constructs from C and C++, Java has different semantics, focusing more on code reuse. These days, Java is owned by and is developed by Oracle Corporation.

Java owes much of its popularity by being there when the Internet started to gain popularity. It was the only language that allowed some sort of dynamic behaviour on browsers through *applets*, small applications that run on browser plugins. The language, however, proved to be useful to much more than that.

The language was revolutionary: it was fast and efficient, it was a lot easier to develop than C or C++ and had a lot of features that became usual in modern languages, such as garbage collection and a rigid structure for modules. However, Java still has some problems: it is verbose, imposes practices that make development of small programs very difficult and lacks some useful abstractions. Anyway, it is a very good language, and its continued success is proof of that. You can find more information on its [website](https://java.com/).


## Comments on the Language
Java mandates that every aspect of a program is structured inside classes, thus making everything an object. While that makes sense for a lot of programming tasks, some procedures end up looking kind of strange, such as when the same operation must be done on different objects. The Fourier Transform is such a case: you can take the Fourier Transform of a lot of different kinds of sequences, but you must encapsulate it in a class.

The language has some constructs to simplify this, but they still look strange. In this case, it is usual to create a class called something like `FastFourierTransformer` with static methods to perform the computation, and pass data to it when needed. While this kind of makes sense, there should be better ways to do it.

The approach I followed in this code is different: since the objective of these programs is to assess the performance of the language, all the code was put in a single class, and the Fourier Transform behaves pretty much as an isolated function. This hurts code reuse, of course, but is the best way to do this in this case (and, of course, if you need to use the functions, you can just copy the relevant parts).

Java's design has an immense concern for program consistency, and do that by disallowing some constructs. There is no multiple inheritance or operator overloading, for example. The reason is that these abstractions can hurt the program: they tend to be overused and can hurt performance. But, even though that might be true for some situations, both are extremelly useful in certain cases. In these programs, for example, operator overloading could allow me to write

```
X[k] = X[k] + Wkn * x[n];
```

which is a lot simpler that what I actually had to write:

```
X[k] = X[k].add(Wkn.mul(x[n]));
```

which looks weird and difficult to read.

This implementation was refactored (check the repository history) to make the complex number and testing methods modules that can be imported. This is, of course, the _right_ way to do things, as code reuse is easier to mantain. It is easy to create Java modules and packages, although there are some details that must be followed that can be a little confusing if you're not used to them. But once you pick up, the worst decision you have to make is how to distribute your classes among files.

## The Code
Given the way that Java works with modules, there are some files and a folder in this directory. Let's start with the folder.

The `fft` folder contains a package with the implementation of the Fast Fourier Transfor methods, the complex numbers mini-library, and the methods to test and measure execution time. These are the files:

1. `Test.java`: this file contains methods to inspect the results of a FFT functions.

2. `Complex.java`: contains the implementation of a `Complex` class that has functions and methods to deal with complex number operations.

3. `FFT.java`: this file contains the implementation of the Fourier Transforms themselves. There are four implementations: `directFT` is the O(N^2) complexity implementation using the definition, it is used for comparison purposes; `recursiveFFT` is a recursive implementation of the FFT for vectors of power-of-2 length; `iterativeFFT` is an iterative implementation of the same algorithm, that preserves memory and execution time by avoiding recursive calls; finally, `recursiveNFFT` is an implementation of the FFT decomposition for composite numbers (that is, the length of the vector is a composite number.

The programs in the root directory are as follows:

1. `MainFFT.java`: main program to run power-of-2 length FFTs a number of times and print the average running time for all implementations.

2. `MainAnyFFT.java`: main program to run composite-number length FFTs a number of times and print the average running time for all implementations.

3. `MainTest.java`: this is a simple program to compute the FFT for vectors of some lengths, using all implementations available on the libraries, and printing the results to check if the results are right. A more correct way to implement this was to actually build a test suite, but that could get too complicated, and I don't think it's needed for a small and simple piece of code as this.

## Compiling and Running
If you use Linux, Java is certainly installed on your system, though you might need to install the development package. On Windows or MacOS, you will probably have to download the packages, but installation is easy. There are also a lot of other development kits, such as GNU Java and others. Feel free to install and experiment with all of them.

I used the OpenJDK runtime environment and compiler. To compile the program, just type:

```
$ javac -cp . MainFFT.java
```

This will generate some `.class` files in the same folder; one of them will be called `MainFFT.class` -- this is the file that should be run. You can do that by typing:

```
$ java MainFFT
```

To compile and run `MainAnyFFT.java`, follow the same steps, just change the name of the file you're compiling and the class you're running. The test program can be built by substituting `MainTest.java` instead.

If you want to, you can use these implementations in your code. However, while these methods are working correctly as far as I was able to test them, it is strongly recommended that you don't, as the purpose of this code is educational only. Instead, I advise you to search for a mature and tested library, such as [FFTW](http://www.fftw.org/), which will certainly give you better results.
