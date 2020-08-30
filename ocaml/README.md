# OCaml Version

This folder contains the OCaml version of the Fourier Transform. OCaml is a functional language designed by Xavier Leroy in 1996 by adding object orientation to Caml, itself a dialect of the older  language ML. The original language was a very simple purely functional language that was very expressive and easy to read. OCaml typing is static, but most, if not all, type declarations can be omited because the compiler is able to infer data type from the operations applied on them. While this is a nice feature, OCaml implementation causes some problems (see below).

Although it is a functional language, OCaml has constructs to deal with imperative commands, such as `for` and `while` loops. It is very possible to write (almost) purely imperative programs using that. While that can be a good thing for people learning the language or developing in a hurry, this can also be a problem, since it is easier to program in OCaml using the functional paradigm.

There are some software out there written in OCaml, but it find its main uses in academy -- as is common with functional languages, since they provide a nice representation for theoretical concepts. A language called *F#*, derived from OCaml, is part of the .NET framework and -- maybe ironically -- the [FFTW](http://www.fftw.org/) library for computation of the Fast Fourier Transform, has a lot of parts written in it. You can find more on the language's [website](http://ocaml.org/)


## Comments on the Language

The main issue that a programmer will find with OCaml is learning to think *functionally*, that is, writing their code in functional fashion. This is a problem with any functional language, of course, but OCaml has a small additional problem: since the language has some support for traditional (*ie*, imperative) programming, most programmers will *default* to the paradigm. The problem is that this can be somewhat convoluted in OCaml, and programming functionally is way easier, once you grasp the concepts.

To me, there was a worse problem: since the compiler needs information to do type inference, OCaml does not allow *any* type mixing -- even the most common operators can't deal with different data types. This means that you can't use a `+` sign to add two float numbers, there is a special operator for that, `+.`. Unfortunatelly, OCaml also doesn't perform *upcasting* of numeric types, because, according to the creators, casting is expensive and error prone. That might be true of *downcasting* (*eg*, automatically converting floats to integers), but I still have to find a situation where upcasting causes any trouble.

This leads to very complicated expressions, which, in my opinion, doesn't help with the language expressivity. Since there are a lot of type mixing in the computation of the Fourier transform (integers, floats and complex numbers), it is common to write an expression like

```ocaml
let w = exp { re=0.; im= -2. *. pi *. (float k) /. (float nx) } in
```

when in other languages (Python, for example) I could write

```python
w = exp(-2j*pi*k/nx)
```

which makes *a lot* more sense. Also, notice that at some points, you *have* to put spaces between operators, or the compiler gets confused. For example, in the code line above, there *must* be a space between `im=` and `-2.`. And, as a last remark, float numbers *must* be written with a trailing decimal dot, or the compiler will think that it is an integer number and refuse to compile.

Things get worse when you have to operate over lists or vectors. This harms code legibility, and since computer programs most of time deal with mathematical computations, legibility of mathematical expressions *should* be a design goal.


## The Programs

There are four programs in this folder:

1. `fft_list.ml`: this file implements three functions: `direct_ft` is straight definition of the Fourier transform, `tr_direct_ft` is the tail-recursive version of the same algorithm, and `recursive_fft` is the recursive Fast Fourier Transform. All of these functions operate over lists, and there is no iterative version, since it would be too cumbersome to do that. These functions deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `fft_array.ml`: in this file, you will find three functions: `direct_ft`, which computes the Fourier transform by definition; `recursive_fft`, the recursive implementation of the fast algorithm, and `iterative_fft`, the iterative in-place fast algorithm. These operate over OCaml arrays, also for power of 2 length;

3. `anyfft_list.ml`: implements two functions: `direct_ft`, which computes the transform by definition, and `recursive_fft`, the recursive fast algorithm, both operating on lists. These functions operate on vectors of any composite length. If the length of the vector is a prime number, it falls back to `direct_ft`, and shows no gain in efficiency;

4. `anyfft_array.ml`: implements two functions: `direct_ft` and `recursive_fft`, with the same sematics as the file above, but operating over OCaml arrays.

OCaml has a `Complex` module to implement operations with complex numbers; this module was used.


## Running

OCaml is probably present and installed in your system if you have a Linux box; if it isn't, it certainly can be found on the package manager. Windows and MacOS users will have to download an executable from the [website](http://ocaml.org/). Installation is easy, and, since there is a package manager in the language toolchain, installation of other modules will be easy too.

OCaml is a compiled language, but programs can be run with a single command (that will compile and run the code). This can be done with:

```
$ ocaml fft_list.ml
```

If you want to generate an executable, just type:

```
$ ocamlc -o fft_list fft_list.ml
```

This will generate an executable file named `fft_list` in the same folder, that can be run with the command:

```
$ ./fft_list
```

You can change `fft_list` for any of the programs above to run them.
