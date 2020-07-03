# C# Version

C# is an object oriented language created by Microsoft in the year 2000. It was developed to be the main language using the Common Language Runtime (CLR), a project to try to unify the libraries using in development of Windows applications. It is object-oriented, general purpose, strong typed and offers a lot of resources common in modern languages, such as automatic memory management. C# was made tightly integrated with *Visual Studio*, a suite of compilers and development tools created by Microsoft. It is possible to write programs that don't use it, but the application was created to give support to the development of applications using the language.

Some people say that C# is an imitation of Java, and they might not be completely off, since the two languages share a lot (and I mean a *lot*!) of characteristics, and their syntax is nearly identical -- and since Java was created about 7 years before, it is difficult to deny. Sharing a lot of characteristics with Java, C# is a language that is easy to write in -- programs written in one language can be converted to the other by changing function names and this or that detail on the syntax -- given that no special resources in each language was used (that was the case with these implementations of the Fast Fourier Transform).

If you are working on a Windows machine, you probably need to instal Visual Studio to compile these programs. If you have a Linux or Unix system, however, you can install the Mono Suite, that implements a lot of the components of CLR, and has a C# compiler.


## Comments on the Language

There isn't much to write about the language: it is what you expect from a Java-derived language. It offers the same structures and has a very similar execution model, which leads to very similar code. In fact, although a `diff` of the two versions shows difference on most lines, those are most comestic: name of functions or modules. C#, however, has operator overloading. Although it is a capacity often shunned by developers, since it can be abused (and not rarelly *is*), it is a nice resource when you are implement algebraic types, such as complex numbers.

That being said, there was no difficulty in building the code, and it works as expected.


## The Programs

There are two programs in this folder:

1. `fft.cs`: this implements `direct_ft`, `recursive_fft` and `iterative_fft`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.cs`: this implements `direct_ft` and `recursive_fft` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. The language already has a module of complex numbers, but implementing them was an experience I thougt it was worth it. It was very easy to write it, and -- with operator overloading -- it turns out very easy to use.


## Compiling and Runing

If you are on a Windows system, you will have to install [Visual Studio](https://visualstudio.microsoft.com/) to compile and run the programs. There is a free version with limited capacity, but it is enough to run these programs. Please, go to the site, download it and follow installation instructions. You will probably have to build a project to compile, or use the command line tools (plase, see the documentation on how to do this). In Linux, you must install the Mono Suite -- the package manager of your system probably has a package for it.

Once the compiler is installed, all you have to do is run the compiler to create an executable. In Linux, this is done in the command line by issuing the command:

```
$ mcs fft.cs
```

This will generate an executable file `fft.exe` in the same folder. To run it, first check that it has the execution permission, and type the command:

```
$ ./fft.exe
```

This will run the functions some times, average the execution time and show a table for comparison. The same commands can be done to compile `anyfft.cs`, just substitute `anyfft` for `fft` in the commands above.
