# Ruby Version

This folder contains the Ruby version of the Discrete Fourier Transform. Ruby is another scripting language that bears a lot of resemblance with Python and others. It was created in 1995 because of the creator's dissatisfation with the way that object orientation was implemented by Perl. The resulting language is a little less free form than Perl, more readable and implements a lot of abstractions that are common in object oriented languages. Particularly, there are a lot of syntax that resembles Smalltalk.

Ruby was not very known, but got a lot of tracking with the creation of Ruby on Rails, a web development framework that turned creation of complex pages really easy. It caused a vast migration from Java to the language, and it didn't take long to receive contributions that made it a good choice for every kind of application. You can find more about the language on its [website](http://ruby-lang.org/)


## Comments on the Language

Script languages are, in general, easy to learn and to write in, but have decreased performance. Ruby is no exception, but it is a lot faster now than it was on my first experiments with it. Recent versions of the language implemented a just-in-time compiler, and that might be the reason its efficiency are improved.

The syntax is very clean and objective, and the resulting code is very easy to read. Ruby designed followed what was called the *principle of least surprise*, in general associated with user interface. It means that, if you know Ruby well (or maybe even not so well), then a line of code should do what you *expect* it to do. Of course, that is an unrealistic expectation, but Ruby does it well.


## The Programs

There are a lot of files in the folder, but two programs are to be run:

1. `fft.rb`: this program implements `direct_ft`, `recursive_fft` and `iterative_fft` functions that compute the Fourier Transform, and runs them a number of times and compares the time spent running the transforms. The functions here can deal only when the vectors to be transformed are of power of 2 length (that is, 2, 4, 8, 16, 32, 64, etc.);

2. `anyfft.rb`: this program implements `direct_ft` and `recursive_fft`, that compute the Cooley-Tukey decomposition algorithm for vectors of composite length (that is, the length is a composite number). If the length of the vector is a prime number, it falls back to `direct_ft`, and shows no gain in efficiency at all.

Although complex numbers are not native to Ruby, there is a very easy to use module that implements complex arithmetic. That module was used.


## Running

Ruby is an interpreted language, and is most likely already installed on your system if you use Linux. If not, you can certainly install it from your package manager. Windows and MacOS users can download the files from the site. To run the program, just type:

```
$ ruby fft.rb
```

To run `anyfft.rb`, just substitute the filename on the command line.
