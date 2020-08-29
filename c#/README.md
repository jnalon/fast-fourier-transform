# C# Version

C# is an object oriented language created by [Microsoft](http://microsoft.com) in year 2000. It was developed to be the main language using the [Common Language Runtime (CLR)](https://docs.microsoft.com/en-us/dotnet/standard/clr), a project that unified the libraries used in development of Windows applications. It is object-oriented, general purpose, strong typed and offers a lot of resources common in modern languages, such as automatic memory management. C# was made tightly integrated with [Visual Studio](https://visualstudio.microsoft.com/), a suite of compilers and development tools created by Microsoft. It is possible to write programs without using it, but, since it was created to give support to the development of applications using the language, bigger projects will benefit from using it.

Some people say that C# is an rip-off of [Java](https://java.com) and was created by Microsoft because they couldn't get control of Java. They might not be completely off, as the two languages share a lot (and I mean a *lot*!) of characteristics, and their syntax is nearly identical -- and since Java was created about 7 years before, it is difficult to deny. But, because of that, C# is a language that is easy to learn and to write in - programs written in one language can be converted to the other by changing function names and this or that detail on the syntax - given that no special resources in each language was used (that was the case with these implementations of the Fast Fourier Transform).

Among the differences between both languages, C# has support for things like functional programming, and operator overloading. It is also easier in C# to deal with functions and methods as first-class objects. On the other hand, Java has a bigger community that develops modules for almost any task, which is not so true for C#.



## Comments on the Language

There isn't much to write about the language: it is what you expect from a Java-derived language. It offers the same structures and has a very similar execution model, which leads to very similar code. In fact, although a `diff` of the two versions shows difference on most lines, those are most comestic: name of functions or modules. C#, however, has operator overloading. Although it is a capacity often shunned by developers, since it can be abused (and not rarelly *is*), it is a nice resource when you are implement algebraic types, such as complex numbers.

That being said, there was no difficulty in building the code, and it works as expected.


## The Programs

There are two programs in this folder:

1. `fft.cs`: this implements `DirectFT`, `RecursiveFFT` and `IterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.cs`: this implements `DirectFT` and `RecursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `DirectFT`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. The language already has a module of complex numbers, but implementing them was an experience I thougt it was worth it. It was very easy to write it, and -- with operator overloading -- it turns out very easy to use.


## Compiling and Runing

If you are working on a Windows machine, you probably need to instal Visual Studio to compile these programs; follow the link above, download and install it. Visual Studio is also available for MacOS, and .NET core is also available for Linux ([check here](https://docs.microsoft.com/en-us/dotnet/core/install/linux)). On Linux and Unix systems, you can also use the [Mono Runtime](https://www.mono-project.com/), that implements a lot of the components of CLR, and has a C# compiler. It probably is available on the package manager for your distribution.

I used Mono on my Linux system. Since C# is a compiled language, you have first to pass it to the compiler to generate an executable. This is done by issuing the command:

```
$ mcs fft.cs
```

This will generate an executable file named `fft.exe` in the same folder. To run it, first check that it has the execution permission, and type the command:

```
$ ./fft.exe
```

This will run the functions some times, average the execution time and show a table for comparison. The same commands can be done to compile `anyfft.cs`, just substitute `anyfft` for `fft` in the commands above.
