# Groovy Version

This folder contains the Groovy version of the Discrete Fourier Transform. Groovy was created by the Apache Foundation in 2007. It is one of the many languages that uses the Java Virtual Machine, but implementing a different set of semantics. The goal of the design was to create a language simple enough to be used as a script language, but powerful enough to be used as a programming language. It is inspired by Python and Smalltalk.

It is in fact very simple, looking, at first contact, as a version of Java without all the verbosity and complexity. It has, however, some features that Java doesn't, such as operator overloading (at this point you probably guessed that I like this), dynamic typing and some data structures and operators. That allows for some idioms that are not seen in its mother language. It also provides support for functional programming, such as currying. You can find more about Groovy in its [website](http://groovy-lang.org/).


## Comments on the Language

Since it drives many inspirations from languages such as Python, Groovy is fairly easy to learn, especially if you already program in Java. The use of the Java Virtual Machine has the benefit of making available the huge and very well tested Java Library; on the other hand, that can be a problem, because you can feel comfortable writing 'Java in Groovy' and getting results, and not using the capabilities of the language itself. While I don't think that this will necessarily lead to worse efficiency (*ie*, slower and bigger programs), that will probably make the development process slower and less efficient. I probably made a lot of these mistakes in these files, as a seasoned Groovy programmer will probably notice.

I actually didn't see any problems with the language, besides (maybe) performance, that is not the best for a JVM language. You will notice, however, that it has a Just-In-Time compiler that gets more efficient during the execution of the program. This might give you some strange results, such as that smaller vectors generate bigger execution times.


## The Programs

There are two programs in this folder:

1. `fft.groovy`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.groovy`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

Besides that, I wrote a small complex library. It is very likely that Java's libraries could be used, but writing this was easier than searching for the docs.


## Running

If you use Linux, Groovy is probably on you package manager, install it from there. If you're using Windows or MacOS, you will have to download, install and configure the packages.

Groovy is a compiled language, but you can run it with the following command:

```
$ groovy fft.groovy
```

You might get a warning before the program runs. This is because a small incompatibility of Groovy with some characteristics of the latest Java Virtual Machine. This doesn't affect the execution, and will probably disappear in the future. To run `anyfft.groovy`, just change `fft.groovy` where needed.
