# Erlang Version

This folder contains the [Erlang](https://erlang.org/) implementation of the Discrete Fourier Transform. The language was created in 1986 by Joe Armstrong, Robert Virding and Mike Williams, researchers working on Ericsson, and it was mainly developed to solve problems in communications. The main problem that the language tried to address was that, in telecommunications and especially in telephony, many processes must happen simultaneously, and thus an easy implementation of isolated and distributed processes are needed.

There is another characteristic of the language that needs to be mentioned: it follows the functional paradigm, which is very uncommon in engineering-aimed languages. The traditional paradigm of programming languages is *procedural*, where each function is a sequence of commands that tells the computer *how* to perform the computation. In functional programming, each function is a sequece of declarations that tells the computer *what* it means.

Let's give an example: if you want to compute the factorial of a number in a procedural language, you create an accumulator and enter a loop where products are computed. The loop is executed and, at the end, the accumulator has the result. In pseudocode:

```
acc = 1
for i = 2 to n
    acc = acc * i
```

In functional programming, you tell the computer that the factorial of a number is that number multiplied by the factorial of its predecessor; you also instruct the computer that the factorial of 1 is 1. Something like this:

```
factorial(1) = 1
factorial(n) = n * factorial(n-1)
```

Of course, most languages *do* provide constructs like this, but they are always a set of instructions. More complex programs (such as the Fast Fourier Transform) makes the distintion very clear.

The multiprocessing aspect of the language is also interesting, since it shows a nice forethinking from the creators. Although multiprocessing already existed when the language was created, it took some time to reach the state we are today: nowadays, it seems that almost every program or application uses multiprocessing or distributed processing at one point. Erlang was developed with this possibility in mind. You can find more information about the language in its [website](https://erlang.org/).


## Comments on the Language

Development in functional languages is hard for two reasons: first, it's a different paradigm, a different way of thinking about the program. Second, creators of functional programming languages are *purists*, for whom the theoretical purity of the paradigm is more important than doing something useful with it. So, using integer numbers along with floating point is *forbidden*, because they're different kinds of objects, and you *must* explicitly cast. Input and output is *forbidden*, because a function that receives the same arguments *must* give the same result (which doesn't happen if you're receiving input from the user of using random numbers, for example). To allow for these operations while still maintaining purity, functional language developers usually devise a series of complicated solutions.

Luckily, this is not true with Erlang. While it follows closely the functional paradigm, it also allows for easy sequencing of steps, if you need them, and input and output of data. While this might make the language *"impure"*, it is a practical solution, and most of the time you want a result, not a theoreticall proof of a concept. It turns out that it is easy to write functional programs with Erlang, but using *"impure"* operations when needed.

I can think of some changes that could be made to the language, however, that might improve it. First, the passing of functions as first class objects. Although this should be a given in a functional language, the process is a little more complicated in Erlang than it should be. You must define an *anonymous* function to pass it as argument to another function, even if the only thing this anonymous function does is call the *named* one (somewhat related to this, I would love to see currying too). Second, I would like to see function nesting, so it would be easier to put loops near where they belong, as this could make explicit code hierarchy and help prevent bugs. As it is now, you must define another recursive function to make a loop, and it cannot be inside the calling function. Third, I would improve on the syntax for list comprehensions and operator overloading. The language also mandates that variable names start with uppercase letters, and functions with lower case letters. It is just a case of personal taste, but I like some freedom at these points.

(Also, Erlang didn't let me use math functions in guards, and that affected at least one function. It was not a problem in this case because it was easy to work around, but that might lead to cumbersome expressions in more complicated cases).

On the other hand, Erlang provides a lot of capabilities and I must confess that I'm surprised. As you can see if you inspect the code, you will see that there was a lot of different ways to implement the same behaviour. While that might seem confusing and some experience with the language is needed to understand the best (and most efficient) way to do things, it makes the language very flexible. Programming in Erlang is easy, and I would love to make something bigger with it someday.


## The Programs

There are two programs in this folder:

1. `fft.erl`: this implements `direct_ft`, whici is the implementation of the algorithm from the definition using list comprehensions; `list_direct_ft`, the same algorithm but using faster list functions; `tr_direct_ft`, the same algorithm with list comprehensions once again, but tail-recursive, which is faster and saves resources; `list_tr_direct_ft`, the tail recursive implementation, but with list functions; `lc_dft`, a one-liner list comprehension version, `recursive_fft`, which is a decimation-in-time algorithm using list comprehensions to split the vector in even and odd samples, and `split_recursive_fft`, which is the same algorithm, but using a `split` external function the samples. There is no iterative version, because iterative commands such as `for` and `while` are not present in functional languages. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.erl`: this implements `direct_ft` from the definition with list functions; `recursive_fft`, a decimation-in-time using list comprehensions to split the samples; and `split_recursive_fft`, the same algorithm, but using an external function to split the samples. The implement the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to the `direct_ft`, and shows no gain in efficiency at all.

Besides the transform functions, both files also implement a small library to deal with complex numbers. Since Erlang cannot overload operators or create new operators, operations with complex numbers are cumbersome. Also, because of the problems with anonymous functions, I had to create functions that return anonymous functions to pass as arguments - if it seems convoluted in the explanation, it's because it is like this in the code.

It is impossible to implement a non-recursive version of the fast Fourier transform in Erlang, because, as a functional language, it doesn't have the commands and constructs to implement loops or to set elements on a list. It *is* possible to emulate the behaviour using recursive functions, but the final code would look to much as a weird and inefficient version of the recursive implementation, so it is not done here.


## Compiling and Running

As far as I know, there is only one implementation of the language, the official one, that can be found on the [homepage of the language](https://www.erlang.org/). It is very likely that the package manager of your Linux distribution, if you are in one, has the packages to install the language. If you are on Windows or MacOS, you will have to download the packages and install. The process is easy and well explained on the page.

Erlang compiles to bytecode, which is then interpreted on its own virtual machine (if you're used to Java derived languages, it is **not** the Java virtual machine, but Erlang's own). There are basically two ways to compile and run the programs. If you want to compile on the command line, just use the command below:

```
$ erlc fft.erl
```

This will generate a bytecode compiled file in the same folder called `fft.beam`. It might complain of unused functions: those are functions that were created to test and verify results, you can safely ignore it. To run the program, issue the command:

```
$ erl -noshell -s fft start -s init stop
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

Don't forget to add the period (`.`) to the end of the line. Again, it might complain of unused functions, just igonre it. You can run the program and see the results by typing:

```
2> fft:start().
```

Don't forget the period. Change `fft` to `anyfft` to see the second program in action.
