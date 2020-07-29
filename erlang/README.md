# Erlang Version

This folder contains the Erlang implementation of the Discrete Fourier Transform. The language was created in 1986 by researchers working on Ericsson, and it was mainly developed to solve problems in communications. There are some highlights on the language, however, that deserve mention: it was created to be distributed -- that is, that programs could run easily in multiprocessed or distributed systems -- and it's functional -- a different paradigm of computation.

The functional paradigm is an interesting one. Programs in most of languages are created as a sequence of steps that must be followed in order to reach a result, and programs are //algorithm// implementations. In functional languages, computations are given by a set of relations between variables and a set of initial conditions. It's very different from what most people think as computation, and it makes a little harder to program in.

The distributed part is also interesting, since it shows a nice forethinking from the creators of the language. Although multiprocessing already existed when the language was created, it took some time to reach the state we are today: nowadays, it seems that almost every program or application uses multiprocessing or distributed processing at one point. Erlang was developed with this possibility in mind. You can find more information about the language in its [website](https://erlang.org/)


## Comments on the Language

Development in functional languages is hard for two reasons: first, it's a different paradigm, a different way of thinking about the program. Second, creators of functional programming languages are *purists*, for whom the theoretical purity of the paradigm is more important than doing something useful with it. So, using integer numbers along with floating point is *forbidden*, input and output is *forbidden*, because that is likely against some law of the functional programming police. However, you *need* to mix integers and reals, and you *need* input and output, and functional languages make everything they can to make it as hard as possible.

Luckily, this is not true with Erlang. Erlang follows closely the functional paradigm, but allows for easy sequencing of steps, if you need them, and input and output of data. It is not perfect, but, at least in those and some other points, the difficulties you have are with the programming part, not struggling against the language to print two lines of text.

I can think of some changes that could be made to the language, however, that might improve it. First, the passing of functions as first class objects; in Erlang, you must define an anonymous function to pass it as argument to another function, even if the only thing this anonymous function does is call the first one (somewhat related to this, I would love to see currying too). Second, I would like to see function nesting, so it would be easier to put loops near where they belong; as it is, you must define another recursive function to make a loop, and it cannot be inside the calling function. Third, I would improve on the syntax for list comprehensions. And, lastly, I would improve on operator overloading.

But, that being said, programming in Erlang is easy. I would love to make something bigger with it someday.


## The Programs

There are two programs in this folder:

1. `fft.erl`: this implements `direct_ft`, whici is the implementation of the algorithm from the definition using list comprehensions; `list_direct_ft`, the same algorithm but using faster list functions; `tr_direct_ft`, the same algorithm with list comprehensions once again, but tail-recursive, which is faster and saves resources; `list_tr_direct_ft`, the tail recursive implementation, but with list functions; `lc_dft`, a one-liner list comprehension version, `recursive_fft`, which is a decimation-in-time algorithm using list comprehensions to split the vector in even and odd samples, and `split_recursive_fft`, which is the same algorithm, but using a `split` external function the samples. There is no iterative version, because iterative commands such as `for` and `while` are not present in functional languages. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.erl`: this implements `direct_ft` from the definition with list functions; `recursive_fft`, a decimation-in-time using list comprehensions to split the samples; and `split_recursive_fft`, the same algorithm, but using an external function to split the samples. The implement the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. Since Erlang cannot overload operators or create new operators, operations with complex numbers are cumbersome. Also, because of the problems with anonymous functions, I created functions that return anonymous functions to pass as functions arguments (got lost? I kind of did).


## Compiling and Running

It is most likely that Erlang is not installed in your system, since it's not a mainstream language. But Linux systems usually have it on their package management system, so install from there. The [homepage of the language](https://www.erlang.org/) has download links for most operating systems.

Compiling and running these programs is easy. First, you have to compile it, which can be done with the command:

```
$ erlc fft.erl
```

This will generate a bytecode compiled file in the same folder called `fft.beam`. It might complain of unused functions: those are functions that were created to test and verify results, you can safely ignore it. To run the program, issue the command:

```
$ ./erl -noshell -s fft start -s init stop
```

To compile and run the `anyfft.erl` file, follow the same steps, just change `fft` to `anyfft` in the commands. Once running, the program will repeat the function calls a certain number of times, and show a table comparing the methods.

Alternatively, you can use the command line tool to load, compile and run the program. To do that, invoque the interpreter:

```
$ erl
```

and, from there, load the desired package:

```
1> c(fft).
```

Again, it might complain of unused functions, just igonre it. You can run the program and see the results by typing:

```
2> fft:start().
```

Change `fft` to `anyfft` to see the second program in action.
