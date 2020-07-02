# Python Version

This folder contais the implementation of the Fast Fourier Transform using the [Python Programming Language](https://python.org). Python is an interpreted script language that found great success in scientific environments. There are a lot of reasons for that: Python is easy and clean, it doesn't get in the way of the solution of the problem, is very readable, has an active and passionate community and -- most of all -- have libraries to accomplish whatever tasks you need.

Most of the systems these days come with Python already preinstalled. That is most true in all operating systems based on Unix or Linux, that includes BSDs flavours, probably all Linux distributions and Mac OS X. Windows users can download and install the binaries or search for a Python distribution, such as Anaconda. In any case, it comes with (or it is very easy to install) some fantastic libraries to deal with scientific data. A complete survey of them, however, is beyond the scope of this document.

It is strongly recommended, however, that you install [NumPy](http://numpy.org/), which is a very powerful library to deal with numerical computations. Some of the scripts in this folder use it.


## The Scripts

The scripts in the folder are described below:

1. `fft.py`: This is the base script, and it has a lot of different implementations to explore some  of the capabilities of the language. The functions in this file implement the Cooley-Tukey decimation-in-time algorithm for vectors of power of two length:

   1.1. `direct_ft`: this function is the base implementation. It is a direct translation of the analysis function, and it is not supposed to be very efficient. It is implemented using Python native lists, and it's very readable and easy to understand;

   1.2. `lc_dft`: this is an implementation of the direct FT using Python list comprehensions. Given the clarity and the conciseness of Python, it is a very readable function, but I urge you to explore the source code and see the comments on the function;

   1.3. `array_direct_ft`: this is the same direct FT, but implemented using NumPy arrays;

   1.4. `recursive_fft`: a recursive implementation of a decimation-in-time algorithm for powers of two. If you consult any book on the subject (I, of course, suggest mine, *"Introdução ao Processamento Digital de Sinais", LTC*), you will see that Fourier Transform has a lot of nice redundancies in its definition, which can be used in a very intuitive way with recursive definitions. This implementation uses Python native lists to accomplish the task.

   1.5. `array_recursive_fft`: the same algorithm, but implemented using NumPy arrays;

   1.6. `iterative_fft`: if you studied a little of Computer Science, you will know that, although recursive definitions can be very natural and intuitive, they might not be the most efficient implementation, since the frequent function calls exert an overhead that uses memory and costs time. There is a way, however, that it can be translated to a iterative version, and this is it. This uses Python native lists;

   1.7. `array_iterative_fft`: the same implementation, but with NumPy arrays.

2. `fft_list.py`: contais the `direct_ft`, `lc_dft`, `recursive_fft` and `iterative_fft` defined above, but in a separate file, so that it doesn't need external modules that are not part of the Python standard distribution. With this, you can try other Python interpreters, such as [Pypy](https://pypy.org), [Jython](http://jython.org) or [Cython](https://cython.org).

3. `fft_array.py`: contais the `array_direct_ft`, `array_recursive_fft` and `array_iterative_fft` above, separated in one file.

4. `anyfft.py`: this is the implementation of the Cooly-Tukey algorithm, but for vectors of composite length (that is, which length is a composite number). It implements the following functions:

   4.1. `direct_ft`: this is the same function described above;

   4.2. `array_direct_ft`: this is the same functions described above;

   4.3. `recursive_fft`: this function is a recursive implementation of the FFT but allowing for vectors which length is a composite number. If it is not, it falls back to the direct definition. This uses Python native lists;

   4.4. `array_recursive_fft`: the same algorithm as the previous function, but implemented with NumPy arrays;

   4.5. `vec_recursive_fft`: the same algorithm, but twiddle factors are computed using NumPy vectorization, which is expected to run faster.

5. `anyfft_list.py`: the same functions `direct_ft` and `recursive_fft` as described above, but in a separated file in case you want to try it with differente Python versions.

6. `anyfft_array.py`: the same functions `array_direct_ft`, `array_recursive_fft` and `vec_recursive_fft` as above, but in a separated file.


## Running

Python is an interpreted language, and there is not much needed to run any program. It is recommended that you go to the command line and issue the following command:

```
$ python fft.py
```

Instead of `fft.py`, of course, you can choose any of the other scripts described above. The script will repeat the computation of the Fast Fourier Transform for various sizes of vectors, and report the average time. The `array` versions will also compute using NumPy's internal implementation for comparison purposes.


## Comments on the Language

Python is a very easy language, and it is incredibly easy to develop for it. Most of equations are directly translated to code, and in general, you can have results in a few hours -- if you never used the language before. You can get quicker results if you are already proficient on the language. I have been using Python for more than 20 years now, and it was the base of every simulation and graphic that I run in this time (including the nice pictures for my book). So, for me, there was no problem in converting the FFT algorithm to Python.

That being said, there are some things in Python that I would like to see changed. First, Python is slow, even if you use external libraries (but that can be mitigated with the use of Pypy, which includes a JIT compiler). There are some syntathic sugar that I would like to see: range represented in "slice" notation (`begin:end:step`), partially evaluated functions using the anonymous variable (`f(_, a, b)`), and some others. But, overall, I'm very satisfied with the language.
