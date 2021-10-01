# Python Version

This folder contais the implementation of the Fast Fourier Transform using the [Python Programming Language](https://python.org). Python is an interpreted script language that found great success in scientific environments. There are a lot of reasons for that: Python is easy and clean, it doesn't get in the way of the solution of the problem, is very readable, has an active and passionate community and -- most of all -- have libraries to accomplish whatever tasks you need.

Most of the systems these days come with Python already preinstalled. That is most true in all operating systems based on Unix or Linux, that includes BSDs flavours, probably all Linux distributions and Mac OS X. Windows users can download and install the binaries or search for a Python distribution, such as Anaconda. In any case, it comes with (or it is very easy to install) some fantastic libraries to deal with scientific data. A complete survey of them, however, is beyond the scope of this document.

It is strongly recommended, however, that you install [NumPy](http://numpy.org/), which is a very powerful library to deal with numerical computations. Some of the scripts in this folder use it.


## Comments on the Language

Python is a very easy language, and it is incredibly easy to develop for it. Most of equations are directly translated to code, and in general, you can have results in a few hours -- if you never used the language before. You can get quicker results if you are already proficient on the language. I have been using Python for more than 20 years now, and it was the base of every simulation and graphic that I run in this time (including the nice pictures for my book). So, for me, there was no problem in converting the FFT algorithm to Python.

Python has a lot of stuff that I like and that I wanted to see implemented in other languages. Things like slice semantics (which are wonderful to use with lists and arrays), operator based constructors (like `[]` for lists, and `{}` for sets, significant increase on readability), and list comprehensions (static languages could generate fantastic code for them). I don't know why these characteristics are not more common in other languages (but then again, I'm not a language designer).

That being said, there are some things in Python that I would like to see changed. First, Python is slow, even if you use external libraries (but that can be mitigated with the use of Pypy, which includes a JIT compiler). There are some syntathic sugar that I would like to see: range represented in "slice" notation (`begin:end:step`), partially evaluated functions using the anonymous variable (`f(_, a, b)` definind a new function), and some others. But, overall, I'm very satisfied with the language.


## The Scripts

The scripts in the folder are described below:

1. `fft.py`: This is the base script, and it just calls functions defined in modules for timing and prints the resulting table. This will call the Cooley-Tukey decimation-in-time algorithms for vectors of power of two length.

2. `fft_pure.py`: This is basically the same as the script above, but only imports the `fft_list.py` and `fft_array.py`, which are native Python libraries that are supported by other Python interpreters suchs as [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org), so you can use them to benchmark the implementations.

3. `anyfft.py`: Base script for the Cooley-Tukey decimation-in-time algorithms for vectors of composite length (that is, the length is a composite number).

4. `anyfft_pure.py`: This is basically the same as the script above, but only imports the `fft_list.py` and `fft_array.py`, which are native Python libraries that are supported by other Python interpreters suchs as [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org), so you can use them to benchmark the implementations.

5. `fft_list.py`: Contains implementations of the algorithms using only Python lists, so that it doesn't need external modules that are not part of the Python standard distribution. With this, you can try other Python interpreters, such as [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org).

   5.1. `direct_ft`: this function is the base implementation. It is a direct translation of the analysis function, and it is not supposed to be very efficient. It is implemented using Python native lists, and it's very readable and easy to understand;

   5.2. `lc_dft`: this is an implementation of the direct FT using Python list comprehensions. Given the clarity and the conciseness of Python, it is a very readable function, but I urge you to explore the source code and see the comments on the function;

   5.3. `recursive_fft`: a recursive implementation of a decimation-in-time algorithm for powers of two. If you consult any book on the subject (I, of course, suggest mine, *"Introdução ao Processamento Digital de Sinais", LTC*), you will see that Fourier Transform has a lot of nice redundancies in its definition, which can be used in a very intuitive way with recursive definitions. This implementation uses Python native lists to accomplish the task.

   5.4. `iterative_fft`: if you studied a little of Computer Science, you will know that, although recursive definitions can be very natural and intuitive, they might not be the most efficient implementation, since the frequent function calls exert an overhead that uses memory and costs time. There is a way, however, that it can be translated to a iterative version, and this is it. This uses Python native lists;

   5.5. `recursive_nfft`: an implementation of a recursive decimation-in-time algorithm, but factored to deal with sequences which length is a composite number (that is, a product of primes). It decomposes using the smallest prime and recursively calling the function on each sub-sequence; if the sequence is not decomposable as such, defers it to the `direct_ft` implementation.

6. `fft_array.py`: Contains implementations of the algorithms using only Python arrays from the `array` module, so that it doesn't need external modules that are not part of the Python standard distribution. With this, you can try other Python interpreters, such as [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org).

   6.1. `direct_ft`: this is the same direct FT, but implemented using Python arrays;

   6.2. `recursive_fft`: the same algorithm, but implemented using Python arrays;

   6.3. `iterative_fft`: the same implementation, but with Python arrays.

   6.4. `recursive_nfft`: the same implementation, but with Python arrays.

7. `fft_numpy.py`: Contains implementations of the algorithms using NumPy arrays. Your local installation might not be able to run this file using [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org), but you can try.

   7.1. `direct_ft`: this is the same direct FT, but implemented using NumPy arrays;

   7.2. `simplified_ft`: this implements the direct FT, but using some of NumPy broadcast power to speed up the operations;

   7.3. `matrix_ft`: this implements the direct FT, but constructing the transform kernel and matrix multiplying it by the input sequence;

   7.4. `recursive_fft`: the same algorithm, but implemented using NumPy arrays;

   7.5. `iterative_fft`: the same implementation, but with NumPy arrays.

   7.6. `recursive_nfft`: the same implementation, but with NumPy arrays.


## Running

Python is an interpreted language, and there is not much needed to run any program. There is a number of different ways to run the program, depending on what you want to inspect. If you want to see all of the algorithms running for power-of-two sequences, go to the command line and type:

```
$ python fft.py
```

The script will repeat the computation of the Fast Fourier Transform for various sizes of vectors, and report the average time. If you want to see the performance of the algorithms with composite-length sequences, then type:

```
$ python anyfft.py
```

You can also run the programs with other implementations of Python. Not every one of them can support modules like NumPy, so you will have to deal with the `pure` versions: `fft_pure.py` for power-of-two sequences, and `anyfft_pure.py` for composite-length sequences.
