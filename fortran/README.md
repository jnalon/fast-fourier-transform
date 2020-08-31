# Fortran 77 Version

This folder contains the Fortran 77 version of the Discrete Fourier Transform. Fortran is one of the very first programming languages, created in 1956 by John Backus, to make coding of numerical computations easier. It's name was firstly written *FORTRAN* (in uppercase letters), from *FORmula *TRANSlator*, that is, a language to translate equations to computer code. It passed through a lot of revisions, the 1977 one being one of the most common in existence.

It is a very old and dated language, and it shows. However strange it may seem, one ought to remember that it was created when there was absolutely no theory on programming languages and, in fact, its creator became one of the first people to study them. It was actually created to be easily translated to machine code while still being readable by a knowledgeable person. At that time, it was quite an acomplishment -- and you *can* actually see that in the program in this folder.

The language is full of idiosyncrasies: the first five columns of each line is reserved for line numbers, and the sixth for control (such as identifying the line as a comment or a continuation of the previous line), so the commands must be written from the seventh column forward. But no line of code can be longer than 72 columns, so you actually have little space to write your code. Restrictions like that are because, at that time, programs were written in punched cards, one code line per card, and they couldn't be longer than 80 columns.

Structured commands didn't exist, and you had to deal with them using some form of `GOTO`, which lead to very ugly code, especially from people who were not used to programming. Fortran saw a lot of revisions and, even though it is an old language, it can look a lot more modern these days. The last revision was in 2018, so not long ago, and added a lot of modern resources, such as concurrent programming.


## Comments on the Language

Who is using Fortran 77 these days? Why should I bother making a version using the language? It's true that Fortran 77 is an outdated and most obsolete programming language, but there are some reasons it is here, most of them personal. The first one is that it was the first language that I learned when I went to college. At that point, I already programmed in BASIC (Microsoft version installed on MSX computers) and Pascal (Borland's Turbo Pascal), and I enjoyed learning a new language. Although, at that time (the 90's), computer languages were *very* different from what they are today, Fortran was *already* a dated language, but it was what we had to learn.

The second reason is that the first time I saw Fast Fourier Transform code, it was written in Fortran. It was hideous, it still gives me nighmares (not really, but you get the point). I'm saying things like using a numbered `RETURN` statement to emulate a `CASE` -- which *already exists* in Fortran. I knew I could do better (but, comparing to that code, any person could).

The last reason is that there are some legacy software around there that might benefit from this. I don't think that any old software written in Fortran 77 around the world which needs a Fast Fourier Transform doesn't already have it, but anyways.

If you don't know anything about old languages, this is a great oportunity to see code written in one of them. It feels old, looks weird and there are a lot of restrictions that look artificial today, but were absolutely needed at the time. I have no other comments about the language, except that, as strange as it may seem, it was incredibly easy to write the code.


## The Programs

There is only one program in the folder:

1. `fft.f`: this program implements `DirectFT` and `IterativeFFT` functions that compute the Fourier Transform. These functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.). There is no recursive version of the code because Fortran 77 can't work with recursive functions. For this exact reason, there is no version of the code to deal with prime factors.

Fortran deals natively with complex numbers, and it is extremely easy to deal with them. This was expected, after all, since the language was created to make numerical computations easier to program.


## Running

There are a lot of Fortran implementations out there. I'll be using the GNU version, `gfortran`, but you can find other compilers around the Internet. Windows and MacOS users will need to do that, in general, installation is easy and straightforward. To compile the program, just type:

```
$ gfortran -o fft fft.f
```

This will generate an executable file names `fft` in the same folder, which can be run by issuing the command:

```
$ ./fft
```
