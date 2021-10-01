# Dart Version

This folder contains the Dart version of the Discrete Fourier Transform. It is an object oriented language created with the goal to power mobile and web applications. It was developed by Google and was first released in 2012. Dart is one of the numerous languages to have its syntax inspired by C, but with a different semantics that borrows heavily from Java, JavaScript and other object-oriented languages. Despite the Wikipedia page not mentioning it, I can see a lot of inspiration from Python too, maybe indirectly from Kotlin.

The language itself is very simple, but very powerful. It is statically typed, but types are inferred most of the time; you can declare a variable as dynamic, and it can hold any kind of object (this is useful for lists of objects, for example). Dart can compile for a variety of different platforms: native, mobile or JavaScript. This makes the development of portable applications easier. Dart is the base of Flutter, a framework to develop for mobile and web that, at the time of this writing, is gaining some traction.

You can find more information on its [website](https://dart.dev/)


## Comments on the Language

If you know any C derived language, you can read and understand Dart code. You will probably understand it a little better if you have knowledge of C++, C# or Java, as there is a lot of these languages put into Dart. Besides the usual way of using `if`s, `while`s and defining classes, Dart has its own peculiarities.

For example, Dart doesn't allow for function overloading. This can be a problem if you are trying to write a class, because usually constructors *need* to be overloaded. They solved this problem with *named constructors*, that you can build different constructors with different -- and descriptive -- names.

Despite not allowing function overloading, you *can* overload operators, and that is a plus. However, you can create only one version of the overloaded operator for each class, which means that, if you need to process different types of objects, you need to resort to some black magic -- that I still don't understand exactly how it's done.

There is also a *cascading operator* for objects, represented by two dots (`..`). You use this to select a method from an object but, instead of getting back the returning value of the method, you receive back the instance. This allows the programmer to cascade multiple calls to the same object without programming specific methods for that.

Finally, the language has a lambda syntax that's really intuitive.

There are two things that I didn't like in the language and that could actually prevent me from using it in projects. First, I think that some commands could be made shorter, and that will improve readability (which, I think, wasn't one of the main goals of the designers).

Second, and this is terrible: *Dart doesn't have a built-in way to format output!*. You can interpolate strings very easily, but if you need to format your output, you have to program it yourself. Granted, this can be easily done, but it is error prone and is such a common task that any language should have -- this *should not* be delegated to the programmer.

(I might be wrong in this second thing, I'm still searching, but couldn't find anything about that).


## The Programs

There are two programs in this folder:

1. `fft.dart`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.dart`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

A small library to deal with complex numbers is also on the files. They were very easy to write, and operator overloading allows for very easy coding of the expressions.


## Running

On my system (Ubuntu), Dart is not installed nor in the package manager, and the editor I use (Kate) doesn't support it. To install it is easy, though, you just have to add a repository -- their webpage has more details on that.

Dart can compile to native, mobile and web, but I tried only the native build. You call it just like an interpreter, and it will compile and run the code (but won't generate an executable file, however). Use the following command:

```
$ dart fft.dart
```

To compile `anyfft.dart`, substitute `fft.dart`. This will run the tests a number of times and print time information on the screen.
