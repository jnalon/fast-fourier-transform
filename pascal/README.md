# Pascal Version

This folder contains the Pascal version of the Discrete Fourier Transform. Pascal is an old language, and its history is well known: Niklaus Wirth designed it in 1970 to use as a tool to teach structured program for his students. At that time, structured programming was the newest method to develop good code, and it seems weird these days that it wasn't common practice. Pascal was based on ALGOL and Simula, but was simpler, and was extremelly successfull.

A **big** part of its success was the emergence of an implementation that was easy to use and not expensive: Borland launched its *Turbo Pascal* compiler, which run in a lot of different operating systems and also featured an integrated development environment (though very rudimentar for today), turning development easier and quicker. To understand why this was a big thing, imagine working on an environment that *doesn't* have windows or multitasking. You had to start your text editor, write your code, save it, get back to the operating system, start the compiler, make notes of every error in compilation, go back to your code in the text editor and search for the lines with errors. With Turbo Pascal, you could compile directly from the editor by pressing a function key; when the error was found it took you back directly to where it was in the code.

With the advent of object orientation, Pascal was quickly adapted into *Object Pascal*, that added classes to the language. To this day, Pascal is used a lot in development of small and big applicatons.


## Comments on the Language

It is difficult for me to decide if Pascal is easy to learn and code or is a good programming language overall: I've been using it for more than 30 years now: Turbo Pascal, in it's third version for CPM operating system, was my first contact with the language. I can say, however, that it is extremely easy to express algorithms in it, and since the Fast Fourier Transform is a very well structured algorithm, it was easy to program.

Pascal share a lot of characteristics with C, but at a little higher level. It emphasizes the creation of good structure and placement of functions, procedures and main code, and good programming practices as a whole. These days, there are a lot of implementations of the language, and you can choose your favorite flavour.

That being said, standard Pascal has something that is kind of a disappointment, at least in comparison with modern languages: that arrays don't have information about their length. This is because, at the time of its creation, language designers were more worried about how the memory was used than to make programmer's life *easier* (memory, at that time, was more scarce that programmer's patience). This is not a big deal, of course, but it had some impact on the programs in this folder.


## The Programs

There are a two files in the folder:

1. `fft.pas`: this program implements `DirectFT`, `RecursiveFFT` and `IterativeFFT` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.pas`: this program implements `DirectFT` and `RecursiveFFT`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `DirectFT`, and shows no gain in efficiency at all.

A small set of functions to deal with complex arithmetic is present on each file. Object Pascal allows for operator overloading, and this made the translation of the equations easier.


## Running

There is a great number of Pascal implementations out there, but I suggest you use [FreePascal](http://freepascal.org/). Its a free software implementation, very complete, and mostly compatible with older Object Pascal compilers (namely, Borland's Object Pascal). To compile the program, just type:

```
$ fpc fft.pas
```

This will generate a `fft` executable file in the same folder. The file can be run by typing:

```
$ ./fft
```

FreePascal has an integrated development environment that can be used to compile and run the programs, making it easier. If you are on Linux, you probably can install from the package manager; download from the webstie for Windows and MacOS. To compile and run `anyfft.pas` just substitute in the command line.
