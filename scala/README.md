# Scala Version

This folder contains the Scala version of the Discrete Fourier Transform. It was created by the Apache Foundation in 2004, with the aims of being a general purpose object oriented functional language, simple to write good code on. It runs on the Java Runtime Environment, and as all of this kind of language, it has access to the vast, complete and mature library of Java modules. That is always a plus for any language.

According to its creators, Scala was designed to address criticisms of Java. It implements some abstractions that are not seen in Java (such as operator overloading and first-class functions), but it has some idiossincrasies, such as the existence of companion objects to classes. It also implements functional paradigm. You can find more about the language on its [website](http://scala-lang.org/)


## Comments on the Language

Scala was designed to be easy to write and read, and it does the job. The programs are very concise, although they can get a little verbose especially when using methods derived from the Java library. The functional paradigm implements very welcome constructs, such as for-expressions which work kind of like list comprehensions from Python and Haskell. In fact, my code has a function that implements the Fourier Transform using them.

The problem with languages based on JRE is that you might feel tempted to write Java code using a different syntax. This can certainly happen with Scala, and it is likely that there are points in my programs where I could use different constructs to get better expressivity or efficiency. But that can be seen as an advantage, because Java bytecode runs very fast.


## The Programs

There are two programs in the folder:

1. `fft.scala`: this program implements `directFT`, `forFT`, `recursiveFFT` and `iterativeFFT` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms. `forFT` is a direct implementation using for expressions, resulting in a very small yet very readable expression. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.scala`: this program implements `directFT` and `recursiveFFT`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `directFT`, and shows no gain in efficiency at all.

A small library of complex number arithmetic was added. The existence of operator overloading allows the expressions to remain readable.


## Running

Scala compiles to Java bytecode. This can be done by issuing the command:

```
$ scalac fft.scala
```

This will generate a lot, **a lot** of files in the folder. Those are probably to facilitate the interface with Java, but the amount of files frightened me out. Of all those files, you want to run the one with the definition of the algorithms, which is `fft.class`. To run, type:

```
$ scala fft
```

To run `anyfft.scala`, just substitute the filename on the command line.
