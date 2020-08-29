# Ruby Version

This folder contains the Ruby version of the Discrete Fourier Transform. [Ruby](http://ruby-lang.org/) is a scripting language created by Yukihiro Matsumoto, and bears a lot of resemblance with Python and others. It was created in 1995 because of the creator's dissatisfation with the way that object orientation was implemented by Perl (and any other scripting language at that time). The resulting language is freeform than Perl, more readable and implements a lot of abstractions that are common in object oriented languages. Particularly, there are a lot of syntax that resembles Smalltalk.

Ruby wasn't well known, and far from a mainstream language, until about 2005, when Ruby On Rails was released, a web development framework that turned creation of complex pages really easy. It was revolutionary at that time, incredibly easy to program, and had a lot of ideas that influenced the creation of similar frameworks for other languages. It caused a vast migration from Java, and it didn't take long to receive contributions that made it a good choice for every kind of application. You can find more about the language on its [website](http://ruby-lang.org/).

The main philosophy in Ruby's design is stated as *Principle of Least Astonishment*, meaning that the language should behave in the less confusing way. This is recognized as a goal difficult to achieve, given that programmers come from different backgrounds, and each of them has a set of expectiations of what should be the behaviour of a language construct. According to the language creator, however, this means that, after you know Ruby well (or even not so well), whatever new features you find won't cause any surprise.


## Comments on the Language

My first experiments with Ruby date to the 1.6 or maybe 1.8 era, and I was incredibly disappointed at how slow it was. I'm glad to find out that this is not an issue anymore: Ruby runs fast (as fast as possible for an interpreted scripting language) and get the job done. Also, its standard and contributed libraries have grown significantly, and it is now easy to develop in.

If you know any programming languages, Ruby won't have any surprise for you at first. Every standard construct (from variables to loops) work exactly as expected, so you will be able to produce code -- even good code -- in your first attempts. Of course (as any other language), Ruby offers more than that. Since one of its goals was to develop an *object oriented* scripting language, you will probably find out a lot of power in writing code in that way.

The programs in this folder might not have a lot of object orientation in it, but that's because the Fourier transform works more as an operator than as an object. This means that it will be better implemented as a set of functions instead of an instance of a class. Fortunatelly, Ruby gives full support for first class functions, so that wasn't a problem.


## The Programs

There are two programs in the folder:

1. `fft.rb`: this program implements `direct_ft`, `recursive_fft` and `iterative_fft` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.rb`: this program implements `direct_ft` and `recursive_fft`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `direct_ft`, and shows no gain in efficiency at all.

Although complex numbers are not native to Ruby, there is a very easy to use module that implements complex arithmetic. That module was used.


## Running

Ruby is an interpreted language, and is most likely already installed on your system if you use Linux. If not, you can certainly install it from your package manager. Windows and MacOS users can download the files from the site. To run the program, just type:

```
$ ruby fft.rb
```

To run `anyfft.rb`, just substitute the filename on the command line.
