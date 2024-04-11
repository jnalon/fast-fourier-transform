# C# Version
C# is an object oriented language created by [Microsoft](http://microsoft.com) in year 2000. It was developed to be the main language using the [Common Language Runtime (CLR)](https://docs.microsoft.com/en-us/dotnet/standard/clr), a project that unified the libraries used in development of Windows applications. It is object-oriented, general purpose, strong typed and offers a lot of resources common in modern languages, such as automatic memory management. C# was made tightly integrated with [Visual Studio](https://visualstudio.microsoft.com/), a suite of compilers and development tools created by Microsoft. It is possible to write programs without using it, but, since it was created to give support to the development of applications using the language, bigger projects will benefit from using it.

Some people say that C# is an rip-off of [Java](https://java.com) and was created by Microsoft because they couldn't get control of Java. They might not be completely off, as the two languages share a lot (and I mean a *lot*!) of characteristics, and their syntax is nearly identical -- and since Java was created about 7 years before, it is difficult to deny. But, because of that, C# is a language that is easy to learn and to write in - programs written in one language can be converted to the other by changing function names and this or that detail on the syntax - given that no special resources in each language was used (that was the case with these implementations of the Fast Fourier Transform).

Among the differences between both languages, C# has support for things like functional programming, and operator overloading. It is also easier in C# to deal with functions and methods as first-class objects. On the other hand, Java has a bigger community that develops modules for almost any task, which is not so true for C#.

## Comments on the Language
There isn't much to write about the language: it is what you expect from a Java-derived language. It offers the same structures and has a very similar execution model, which leads to very similar code. In fact, although a `diff` of the two versions shows difference on most lines, those are most comestic: name of functions or modules. C#, however, has operator overloading. Although it is a capacity often shunned by developers, since it can be abused (and not rarelly *is*), it is a nice resource when you are implement algebraic types, such as complex numbers.

That being said, there was no difficulty in building the code, and it works as expected.

Modularizing the program was very easy. It turns out that it suffices to split the functions over appropriated files and let the compiler find the dependencies. You just need to make sure that every module is included in the compilation command.

One of the things that I'm also experimenting with these implementations is the default documentation method. In most implementations, [Doxygen](https://www.doxygen.nl/) comments are enough. However, C# has a standard way to write documentation comments that is -- unfortunately -- based on XML. It is very simple, but it's an uneeded hassle to worry about opening and closing tags (which can also lead to errors if you forget to close a tag and all the problems that arise when using XML). And, considering that developers usually _don't like_ to write documentation anyways, it seems to me that this works against the intent.

## The Code
There are a few files in this folder that implement functionality that is shared among the main programs. The shared modules are the following:

1. `Test.cs`: this file contains methods and functions to inspect the results of a function that implements a Fourier Transform, and also to measure its execution time.

2. `Complex.cs`: contains the implementation of a `Complex` class that has functions and methods to deal with complex numbers operations.

3. `FFT.cs`: this file contains the implementation of the Fourier Transforms themselves. There are four implementations: `DirectFT` is the O(N^2) complexity implementation using the definition, it is used for comparison purposes; `RecursiveFFT` is a recursive implementation of the FFT for vectors of power-of-2 length; `IterativeFFT` is an iterative implementation of the same algorithm, that preserves memory and execution time by avoiding recursive calls; finally, `RecursiveNFFT` is an implementation of the FFT decomposition for composite numbers (that is, the length of the vector is a composite number.

The programs in the root directory are as follows:

1. `MainFFT.cs`: main program to run power-of-2 length FFTs a number of times and print the average running time for all implementations.

2. `MainAnyFFT.cs`: main program to run composite-number length FFTs a number of times and print the average running time for all implementations.

3. `MainTest.cs`: this is a simple program to compute the FFT for vectors of some lengths, using all implementations available on the libraries, and printing the results to check if the results are right. A more correct way to implement this was to actually build a test suite, but that could get too complicated, and I don't think it's needed for a small and simple piece of code as this.

## Compiling and Runing
If you are working on a Windows machine, you will probably need to instal Visual Studio to compile these programs; follow the link above, download and install it. Visual Studio is also available for MacOS, and .NET core is also available for Linux ([check here](https://docs.microsoft.com/en-us/dotnet/core/install/linux)). On Linux and Unix systems, you can also use the [Mono Runtime](https://www.mono-project.com/), that implements a lot of the components of CLR, and has a C# compiler. It probably is available on the package manager for your distribution.

I used Mono on my Linux system. Since C# is a compiled language, you have first to pass it to the compiler to generate an executable. This is done by issuing the command:

```
$ mcs MainFFT.cs Complex.cs FFT.cs Test.cs
```

This will generate an executable file named `MainFFT.exe` in the same folder. To run it, first check that it has the execution permission, and type the command:

```
$ ./MainFFT.exe
```

This will run the functions some times, average the execution time and show a table for comparison. The same commands can be done to compile `MainAnyFFT.cs`, just substitute `MainAnyFFT` for `MainFFT` in the commands above. The same steps can be followed to compile and run `MainTest.cs`.

If you want to, you can use these implementations in your code. However, while these methods are working correctly as far as I was able to test them, it is strongly recommended that you don't, as the purpose of this code is educational only. Instead, I advise you to search for a mature and tested library, such as [FFTW](http://www.fftw.org/), which will certainly give you better results.

