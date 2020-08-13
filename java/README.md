# Java Version

This folder contains the Java version of the Discrete Fourier Transform. Java is an object oriented language created by Sun Microsystems in 1995. It was designed with the goal of allow developing portable applications -- it's *motto* was *\"write once, run everywhere\"*. Although it inherited a lot of syntax constructs from C and C++, Java has different semantics, focusing more on code reuse. These days, Java belongs to and is developed by Oracle Corporation.

Java owes much of its popularity by being there when the Internet started to gain popularity. It was the only language that allowed some sort of dynamic behaviour on browsers through *applets*, small applications that run on browser plugins. The language, however, proved to be useful to much more than that.

The language was revolutionary: it was fast and efficient, it was a lot easier to develop than C or C++ and had a lot of features that became usual in modern languages, such as garbage collection and a rigid structure for modules. However, Java still has some problems: it is verbose, imposes practices that make development of small programs very difficult and lacks some useful abstractions. Anyway, it is a very good language, and its continued success is proof of that. You can find more information on its [website](https://java.com/)


## Comments on the Language

Java mandates that every aspect of a program is structured inside classes, thus making everything an object. While that makes sense for a lot of programming, some procedures end up looking kind of strange, such as when the same operation must be done on different objects. The Fourier Transform is such a case: you can take the Fourier Transform of a lot of different kinds of sequences, but you must encapsulate it in a class.

The language has some constructs to simplify this, but they still look strange. In this case, it is usual to create a class called something like `FastFourierTransformer` with static methods to perform the computation, and pass data to it when needed. While this kind of makes sense, there should be better ways to do it.

The approach I followed in this code is different: since the objective of these programs is to assess the performance of the language, all the code was put in a single class, and the Fourier Transform behaves pretty much as an isolated function. This hurts code reuse, of course, but is the best way to do this in this case (and, of course, if you need to use the functions, you can just copy the relevant parts).


## The Programs

There are two programs in this folder:

1. `fft.java`: this implements `directFT`, `recursiveFFT` and `iterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.java`: this implements `directFT` and `recursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `directFT`, and shows no gain in efficiency at all.

Besides that, I wrote a small complex library. It is very likely that Java's libraries could be used, but writing this was easier than searching for the docs.


## Running

If you use Linux, Java is certainly installed on your system, though you might need to install the development package. On Windows or MacOS, you will probably have to download the packages, but installation is easy. There are also a lot of other development kits, such as GNU Java and others. Feel free to install and experiment with all of them.

I used the OpenJDK runtime environment and compiler. To compile the program, just type:

```
$ javac fft.java
```

This will generate some `.class` files in the same folder; one of them will be called `fft.class` -- this is the file that should be run. You can do that by typing:

```
$ java fft
```

To run `anyfft.java`, just change `fft.java` where needed.
