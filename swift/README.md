# Swift Version
This folder contains the implementation of the Fast Fourier Transform using the Swift Programming language. Swift is a general purpose programming language created by Apple. It was first released in 2014 as a replacement for Objective-C in the development of applications both for iOS and MacOS. It compiles to LLVM, which makes it very portable.

The language has a standard C-based syntax, and it has the standard language constructs, but with a few characteristics that were added to ease development and make programs safer. For example, it has automatic memory management, which _in general_ makes for safer programs, and type inference, which makes the source code less confusing. Also, Swift support _closures_, which are self-contained blocks of code that can be treated like objects.

Since it was designed to take the place of Objective-C in the development of software for Apple products, it inherits a lot of its characteristics and, more importantly, its libraries, _especially_ those ones used to develop for Apple's operating systems. So, I believe, any programmer who is familiar with the development using Objective-C will be at home with Swift.

## Comments on the Language
As I said in the introduction, Swift uses C-like syntax. With just a few exceptions, all modern languages follow the standard, because it makes the learning curve a little gentler. However, this has the side effect of inducing programmers to _"make a C program using Swift"_. By this I mean that your first programs will probably use the simpler constructs (like standard decisions and loops) and won't use the full power of the language. And while the Fast Fourier Transform is a very simple algorithm, I'm probably guilty of doing this on this project, and this impression was stronger when using Swift.

Swift also has one of the features that I like most in modern languages: it is _very_ easy to initialize standard types. For comparison, in C, you have to define the size of the memory used by an array, maybe dynamically allocate it, and if that's the case initialize it by hand. In Python, on the other hand, there are syntax constructs to do that. Swift is more like Python than C in this point, because each standard data structure has a specific syntax and automatically manages the memory.

However, with that being said, Swift has a feature that I _don't_ like in languages, which is the explicit casting of numbers. I understand that bugs relating to miscasting can be very difficult to track and languages should protect against it. However, I also feel that promoting numbers in certain cases won't be a problem. From my experience, miscasting is a problem when you're dealing with data on which bits and their position are relevant (say, when receiving data from an Internet connection), but when you want to take the logarithm of a number, it is _expected_ that your result will be a floating-point number. Casting should be automatic at these cases.

Swift also has a somewhat strange behavior in function definition and calling: parameters to the functions _must_ be explicitly denoted. For example, if I define a function `f` which receive an argument `a`, then in the function call I _must_ identify the parameter: `f(a: 1)`. While this is _possible_ in most languages, and even _desirable_ in many cases, this is the first time I've seen it being mandatory. It is possible, however, to override this by using argument labels -- which are a first for me too.

An _argument label_ is a kind of an alias you can give to an argument in a function declaration. It helps function calling to be made clearer. For example, if you define an argument `n` with a label `count`, you can call the function using the label (_eg._ `f(count: 5)`), which makes it clear what the argument is, but you can use `n` in the internal function definition. That can be good for documentation and keeping the code concise. However, I don't see a lot of gain over just calling your argument `count` to start with. So, I'm not sure if I like this or not.

I also liked the _"anonymous variable"_, represented by an underscore (`_`). By using this variable in certain points, you inform the compiler that you do not care about the result of that operation. It can be used when deconstructing structured data to ignore some fields, to ignore the argument label in function declarations, and to ignore the return value of a function (supressing a warning from the compiler). I used all these three possibilities here, but I'm sure that there's more.

Besides all these characteristics, programming in Swift is standard. I'm pretty sure that there are a lot of other features that make the language more powerful, but I don't think they could be of any use here (except, maybe, for parallel processing capabilities).

## The Code
There are some files in this folder. I'm pretty sure that there is a very easy way to modularize Swift code more than I did here, and I'll sure be looking into it in the future. But, for now, there are only two files:

1. `fft.swift`: this file implements and runs a series of tests of the implementations of the DFT. The implementations are a direct FT, computed by the definition, a recursive version of the FFT and an iterative version of the FFT. In this file, all implementations expect an input vector with a power-of-two size. No checks are made on this, and if the input vector do not follow the restriction, an error will appear.

2. `anyfft.swift`: this file implements and runs a series of tests of the implementation of the DFT for composite numbers, _ie._ vectors which size is a product of prime numbers. These implementations can deal with vectors of any size, but if the size of the input vector is not a composite number, the computation will be deferred to the direct FT, and there will be no efficiency gain.

## Compiling and Running
Since the development of the language is controlled by a single entity, compiling and running the code probably work the same in every environment. Please, be sure that the Swift language is installed (check [their website](http://swift.org/) and follow the installation instructions). I tried this on a Linux box and it worked without problems, so I'll give the instructions to run it from the command line. If you need to incorporate the libraries in your code, please, check your IDE for details.

To run the code, just type the following command:

```
$ swift fft.swift
```

The program will run, execute the tests and return the average execution time. If you want to run the composite-number version, just use the command:

```
$ swift anyfft.swift
```
