# Groovy Version

This folder contains the Groovy version of the Discrete Fourier Transform. It was designed by James Strachan for the Apache Foundation in 2007 as a programming language that could also be used as an scripting language. It runs on the Java Virtual Machine and, because of this, it has access to the vast Java library of modules. It was inspired by features present in dynamic languages as Python and Smalltalk.

As every language based on JRE, Groovy has the disadvantage that programmers might *program Java using Groovy syntax*, but Groovy has this problem taken a little further, because much of Java's syntax is valid Groovy. But Groovy has features that are not present in Java, such as dynamic typing, operator overloading and more native data structures with their own syntax. Groovy also has an interesting *safe operator*, given by `?.`, that checks if an object exists before calling its methods.

Groovy also somewhat supports the functional paradigm. It is not a fully functional programming language, but there are abstractions such as currying, the process of creating functions by partial application, and other featurer that become clear when dealing with collections (lists, maps, sets, and more).

Besides all that, Groovy programs can be executed as if they were scripts. When invoked, Groovy compiles the script to bytecode and runs it on Java's Virtual Machine automatically, but in a process that is transparent to the programmer. You can find more about Groovy in its [website](http://groovy-lang.org/).


## Comments on the Language

Groovy was designed kind of a *Java with less verbosity*, and with a lot of ideas taken from Python and other scripting languages. The main consequence of those choices is that Groovy is fairly easy to learn -- but the dependency on Java, as said above, can lead to write Java equivalent code, but using Groovy syntax. In fact, I'm pretty sure that I made this mistake at some points in these programs.

I liked Groovy's take on closures, though I hadn't used them here. I'm positive, however, that there is a one- or two-liner way to write the Fourier Transform in its direct form using them, with a result much like list comprehensions in Python, Erlang and other languages. I didn't try that yet, but I might get back to it at some point. Closures also are a nice way to pass functions as arguments to functions.

The performance, however, leaves something to be desired. In fact, at some points I got weird results, such as the Fast Transform taking as long as the direct form. I couldn't find any problems with the algorithms, so I suppose that this might be an effect of Groovy's internal: either the function calling process is very slow, or, if it's Just-In-Yime compiler doesn't perform well. In any case, that can be a problem.


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
