# Kotlin Version

This folder contains the Kotlin version of the Discrete Fourier Transform. Kotlin was created in 2011, and is one of the numerous languages that run over the Java Runtime Environment. As such, it benefits from the immense library of modules from, which is always a good thing. Contrary to most but a few of the JRE languages, Kotlin is finding some success because it can be used for development of mobile applications, and also because it can be compiled to Javascript, allowing for easy web development.

Besides allowing access to the extensive library of Java modules, Kotlin tries to make things simpler. Static functions can exist on their own (not needing to be encapsulated by a class), classes can be easily extended without subclassing, there is a lot of type inference in compilation, some objects can be deconstructed, and more. This makes it easier than Java, makes code more expressive, and no power is lost in this.

The downside of being a JRE language, as always, is that people might be tempted to write Java using a different syntax. Kotlin, however, is very simple, and even if you do this, you end up with efficient programs that are fairly readable. Advocates of the language claim that it can be seen as Java if all the verbosity and complexity was removed (Kotlin is not, however, the only JRE-derived language that claims this). You can find more information on its [website](https://kotlinlang.org/)


## Comments on the Language

If you know Java, you probably already know Kotlin -- as said above, it *can* be seen as Java with the complexity stripped down. However, it has different resources that can make development easy, especially on those cases where agile development is needed. Among its features (that I personally like) are first-class functions and operator overloading.

Kotlin *also* allows for functional programming of some kind. It has pattern matching and definition of functions by *what* they do, not by *how* they do. This allows for easy expression of solution to problems and makes it easier to develop creating less bugs. Kotlin, however, avoids the purity *"problem"* by using eager execution and allowing mixing of procedural sequences in the functions. My implementations in Kotlin are *not* functional, but they could be easily transformed into it, given the flexibility of the language.

It needs to be mentioned that languages that mix functional and imperative programming might find that one of the two paradigms, although available, is not used. I think Kotlin is a language like that: being heavily inspired on Java makes it easier to write imperative programs that runs faster, and the functional side ends with little use. It might be, in case of Kotlin, that functional programming *does* give some advantage over imperative, but I didn't see how it could be.

As a final note, the Kotlin compiler is *extremely* slow. The quite small programs in this folder took about eight seconds to compile. While it doesn't seem a lot, I imagine that this could be a problem for larger programs.


## The Programs

There are two programs in this folder:

1. `fft.kt`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.kt`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

A small library to deal with complex numbers is also on the files. They were very easy to write, and operator overloading allows for very easy coding of the expressions.


## Running

On my system (Ubuntu), Kotlin is not installed nor in the package manager. I had to install `snap` and install the language from there. Other Linux-based systems might work in a different way, so check you package manager first. If it isn't there, you can probably install from `snap` (if available), or download from the website. Another possibility is install *Android Studio*, since it comes with a Kotlin distribution (although, maybe, compiling and running these programs might be a little more complicated). As always, Windows and MacOS users will have to download from the site.

Kotlin code is compiled to Java bytecode. The following command does the trick:

```
$ kotlinc fft.kt
```

This will generate some `.class` files (and other files that, probably, are useful for interfacing with Java). Strangely enough, it doesn't generate a `fft.class` to run, but a file with a name that makes it clear that it comes from Kotlin, in this case `FftKt.class`. To run this file, type:

```
$ kotlin FftKt
```

To compile `anyfft.kt`, substitute `fft.kt` in the first command. To run, call the `AnyfftKt.class`.
