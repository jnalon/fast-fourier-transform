# Go Version

This folder contains the Go version of the Discrete Fourier Transform. Go is a general purpose compiled language created by Robert Griesemer, Rob Pike and Ken Thompson at [Google](https://google.com/) in 2009 with the intent to overcome problems that were common in languages such as C and Java, and support features such as multiprocessing. Its syntax is based on C, and it has lots of features of modern languages, such as automatic memory management. You can find more information at their [website](https://golang.org/)

Code hierarchy in Go, however, is very simplified; object orientation, for example, although possible, is not structured with *classes* - you can still program using object orientation philosophy, but constructs are different (based on what the developers call *interface types*, but different from Java's interfaces). Go also doesn't allow for nested functions, which *might* be useful in some cases. These have the following effect: it is extremely easy to produce Go code from the start, but it can take a little while to understand the techniques to produce the *best* Go code. That can be viewed as an advantage, as people learning the language can produce good code from the beginning. As they get experienced with the language, code might get better.

Go also has features to enforce a strong hierarchy of modules and packages. This is certainly good for bigger projects, but can be overkill for simpler projects. But I think the strongest feature is the support for concurrency in the form of coroutines (called *goroutines*) that implements multiprocessing differently from traditional languages; constructs such as locks *are* present to be used when needed).


## Comments on the Language

If you are familiar with procedural languages, it won't be difficult to start producing code in Go. The syntax is similar to C, and the structured commands are basically the same. Bigger projects will, of course, demand that you understand advanced features to be more effective. In these programs, I didn't use anything a lot different from what is done in C or Java, but I suppose that even a simple algorithm such as the Fast Fourier Transform *could* gain a lot from the use of coroutines.

I found array declaration kind of weird, with squared brackets *before* type name, which is unusual and, frankly, looked like Go designers were only trying to be different. The assignment/declaration operator (`:=`) is welcome, since you don't need to clog your code with type declarations. I would like to see initializers for vectors or a better syntax for slices. Also, Go doesn't promote *up* numeric types, even when safe (*ie*, convert integer to float when needed), so you have to add a lot of type casting. That makes the code look ugly at points.

Go doesn't promote *up* number types, and that leads to cumbersome syntax at some points. In the first version of these programs, I would compute a complex exponential, and needed to do something like

```
X[k] = X[k] + x[n] * CExp(-2*pi*math.Pi*float64(k)*float64(n)/float64(N))
```
and that was because I wrote the `CExp` function to make things simpler. You can still find some lines of code such as

```
rn := int(math.Ceil(math.Sqrt(float64(n))))
```
which, frankly, is silly. Promoting up is, as far as I can tell, a safe operation, and there is nothing to be gained by disallowing it -- if you need, cast explicitly *down*.

Besides that, you won't find a lot of surprises, which is probably good: new programmers won't feel *intimidated* by Go. With time, you will learn the best capabilities of the language.


## The Programs

There are two programs in this folder:

1. `fft.go`: this implements `DirectFT`, `RecursiveFFT` and `IterativeFFT`, run them a number of times and compare the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.go`: this implements `DirectFT` and `RecursiveFFT` with the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `DirectFT`, and shows no gain in efficiency at all.

Go has a complex number library. This is fortunate because, since the language doesn't support operator overloading, complex arithmetic code would probably look very cumbersome.


## Running

If you use any flavor of Linux, Go is probably already installed in your system. If it's not, then you certainly can find it in your package manager. If you're using Windows or MacOS, you can [download](https://golang.org/dl/) from their website and install it easily.

Go is a compiled language, and everything related to translation and execution are done the command `go`. To compile and run the program, just type:

```
$ go run fft.go
```

The program will start and run. To run `anyfft.go`, just change `fft.go` to `anyfft.go` where needed.
