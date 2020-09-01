# Frink Version

This folder contains the Frink version of the Discrete Fourier Transform. It was created by Alan Eliasen with the purpose of making physical calculations easy. It's not a mainstream language, and I happened into it when, one day, I was bored and looking for programming resources on an Android tablet. Frink was one of the languages I found that could be used with its own command line and interpreter, making it easy to type quick programs and thus easy to learn. I found out later that Frink can be executed in any environment that has a Java interpreter.

Frink called my attention because of one specific thing: it is the first -- and as far as I know, the only one -- to include manipulation of units of measure embedded in its syntax, and that made me curious. You could write something like:

```
if (1 km == 1000 m) ...
```

and it will automatically make the conversion (this one will result `true`).

In fact, after knowing Frink, I thought how *weird* it is that **no other** language has unit manipulation, since most of them are used and many were created to deal with physics simulations. Maybe it's not a heavily used feature, but it is simple to implement anyway.

Frink is easy to use and learn. It runs on the Java platform, but I wasn't able to find out if it is only an interpreted language //written// in Java or if it compiles to Java bytecode. Anyway, you need to download the files and run, so it doesn't make a strong difference in the end. You can find more about the language on its [website](http://frinklang.org/)


## Comments on the Language

Besides all that, Frink does what other languages do. It runs over the Java platform, it has the usual branching, looping and function commands. It has some requisites (function arguments are specified with brackets instead of parenthesis, block opening brackets should always appear on new lines, and others) but it is incredibly easy to write on. There was absolutely no difficulty in writing the Fourier Transform functions.

That being said, I think it's difficult to me to find out -- at least from the small programs I've written -- if it could be used in larger projects. I think it can't because of some reasons: first, I doesn't seem to generate a compiled `.class` file that could be shared with Java projects, or a machine code module to link with other languages. That means that you'd probably should embed the interpreter with your code, and that seems an unneeded overhead. Frink is nice and has nice capabilities, but that's not something that couldn't be easily programmed in other languages, and including it on projects doesn't seem to bring that much gain.

But, and I emphasize that, you lose nothing learning it. You can do it in a few days, and you might gain some insights. In fact, I think Frink is one of those languages that deserved a little more attention.


## The Programs

There are two programs in the folder:

1. `fft.frink`: this program implements `directFT`, `recursiveFFT` and `iterativeFFT` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms, and also with the internal implementation of the FFT available in the language. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.frink`: this program implements `directFT` and `recursiveFFT`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `directFT`, and shows no gain in efficiency at all.

There is no need to implement a complex library, since complex numbers are transparently handled by the language. In fact, complex numbers are only another type of numeric data, and the language does the conversions when needed.


## Running

To run the programs, you first need to download and install it. The procedure is somewhat convoluted, and not straightforward, but it can be done (of course). As it depends on the Java Runtime, you first need to install it on your system, but that probably isn't an issue, as Java is probably present or easily installed in any or other way. You can find more information on the language's [website](http://frinklang.org/).

I found that it was easier to download the interpreter `.jar` file and run it locally without installing anything. The main page on the website has the link to the file and how to run it. In my system, I did:


```
$ java -cp frink.jar frink.parser.Frink fft.frink
```

There is a `frink` shell script to run on Linux/Unix systems and a `frink.bat` to run on Windows systems, but they didn't work for me. To run `anyfft.frink`, just change the name of the file in the command line.
