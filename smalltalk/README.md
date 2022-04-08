# Smalltalk Version

This folder contains the implementation of the Fast Fourier Transform using the Smalltalk Programming Language. Smalltalk was a revolutionary language created in 1972 by Alan Kay and others. It was the result of the research from Kay in the infamous Xerox Palo Alto Research Center.

It is considered one of the first, if not *the* first, completely object-oriented language -- so, it's needless to say that it influenced every other object oriented language that exist. It was not the first language to have classes or objects, but it is certainly the first one to make these the central point of the design. The goal was to build a language that was closer to how humans think than how computers work. Ohter central point in the design was that the language was part of a complete graphical environment that could be extended and tinkered with.

There are a lot of different implementations of Smalltalk currently. Some are free software and actively developed, such as Gnu Smalltalk, that works kind of a script language, and Squeak, a complete graphical environment. Unfortunatelly, Smalltalk, although praised by its design and elegance, was not a commercial success with the magnitude that C, Java or Python see these days.


## Comments on the Language

The object-oriented feature permeates the whole design of the language. Everything is implemented as an object that responds to messages, and that's basically all there is to it. In fact, even flow control structures common in other languages (such as `if` and `while`) are implemented as messages sent to boolean objects. This makes it quite different from other languages, and a different kind of reasoning.

Because of this, developing the programs was not easy at first. My first difficulty was to understand where to locate the methods that I wanted to write. Since the class hierarchy can be huge, it can be daunting to find what you can change. But, once you understand that, in Smalltalk, you can extend existing classes -- a philosophy of thought that the language design favors -- all was needed was to check the library documentation and choose what class to use. Once this threshold was beaten, programming in Smalltalk was surprisingly similar to programming in Python -- most likely because the object model of Python is based on Smalltalk.

Smalltalk, however, being an interpreted dynamic language, is not fast. In fact, it might be the slowest implementation in this repository. Of course, I'm not a specialist in the language, and certainly there are a lot of ways to improve the code; I don't feel, however, that it will gain a lot of speed.

Overall, I liked working with Smalltalk. It was easy to understand and learn (at the level that I needed to do these programs), and it was fun to be forced to think differently. I didn't like, however, the Gnu Smalltalk version a lot -- I thought that the error messages didn't help a lot, and there are some corner cases in the standard objects that you have to figure out by yourself (for example, using `Float` usually leads to a *division by zero* error, most likely because a `Float` is actually a rational; I had to find out for myself that using `FloatE` solved this problem).


## The Programs

There are two programs in this folder:

1. `fft.st`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.st`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

A small library of complex numbers (based on the example available in the Gnu Smalltalk tutorial) is included in both files, since complex numbers are not standard on the language. It was a nice example of subclassing on the language. The Fourier Transforms are extensions to the `Array` class.


## Running

These programs were written and tested against Gnu Smalltalk, but they are very simple and will probably work on other distributions too. In Gnu Smalltalk, you can run them by issuing the commands:

```
$ gst -f fft.st
```

to run the `fft.st` script, and

```
$ gst -f anyfft.st
```

to `anyfft.st`. The `-f` switch is needed to run the programs in user mode; otherwise, the virtual machine will display debug and warning messages that will mess with the output.
