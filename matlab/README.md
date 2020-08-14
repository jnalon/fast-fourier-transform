# Matlab/Octave Version

This folder contains the Matlab version of the Discrete Fourier Transform. Matlab is an old acquaintance of people who need to deal with scientific and numerical computations for a while. It was created first as an interface to a Fortran package that dealt with linear algebra called LINPACK, but it soon became clear that it had commercial potential, and a new company, called Mathworks, was created. These days, more than a programming language, Matlab became a fully-integrated development environment with lots of modules (called toolboxes) to deal with basically any kind of computational need. You can find more about the language in its [website](http://mathworks.com/). Matlab, the software, also has a graphical programming environment called Simulink.

Given the success of the language, and its somewhat prohibitive price, clones were created. GNU project launched [Octave](http://octave.org), which is mostly compatible with the language and, these days, it also has a complete development environment. [Scilab](http://scilab.org/) is also strongly based on Matlab, and any person who can program in one of these environments are able to program in the other. Both clones are free and open source.


## Comments on the Language

Since its inception, Matlab was created to *deal away* with the complexity: to give access to powerful libraries with easy convertion from equation to code. This can be seen as soon as one starts to learn the language, *especially* if their main goal is to do numerical computation. Syntax is clean, and there are a lot of functions that perform complicated operations, which makes it easy to concentrate in *what* the code should do, not on *how* it shouldd do it.

Of course, it has its problems. For example, each function must be in a separated file if you want to call them externally. This can make development of big programs very difficult, since tracking each file can be cumbersome. Also, it was not created to be *multipurpose*, so you will find it hard to do something different from scientific computing. Object orientation and structures are part of the language, but they were probably done as an afterthought, and using them is not straight forward (although not difficult).


## The Programs

There are a lot of files in the folder, but two programs are to be run:

1. `fft_main.m`: this program calls `direct_ft.m`, `recursive_fft.m` and `iterative_fft.m`. Each of these files has a function (of the same name, as it is mandatory for Matlab) that computes the Fourier Transform according to the respective algorithm. The `fft_main.m` program runs them a number of times and compares the time spent running the transforms, and also the internal Fast Fourier Transform implementations. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft_main.m`: this program calls `direct_ft.m` and `recursive_anyfft.m` and `vec_recursive_anyfft.m`. Each of these files implement a function that runs the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). `vec_recursive_anyfft.m` implements a vectorized version of the recursive Fourier Transform. If the length of the vector is a prime number, it falls back to `direct_ft.m`, and shows no gain in efficiency at all.

Besides that, a file called `time_it.m` implements a function to measure execution time. There was no need to create a complex number implementation, since these are already part of the language.


## Running

You can run these algorithms with any of the programs above; I tested it with Octave and Scilab in Linux, and both run without modification. To install any of these, just look into your package manager; Windows and MacOS users might have to download directly from the site. If you want to use Octave, just issue the following command:

```
$ octave fft_main.m
```

To run `anyfft_main.m`, just change the file name. None of the other files are to be run separately, they're part of these programs.

If you want to run from Matlab, things can get a little complicated: not because it is difficult, but because you will have to acquire a license (which is expensive) or install a trial version (which is limited). But you might have access in an university or company laboratory. In that case, you will need to run the programs from the Matlab command line (this will work with Octave and Scilab command lines also). Start the environment, change to the directory of the programs, and type in the command line:

```
> fft_main
```

and the program will run. Of course, if you type `anyfft_main`, this program will run.
