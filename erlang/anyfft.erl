%% -*- coding: utf-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- Python 3 Version
% This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
%
% José Alexandre Nalon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run, enter the command line interpreter and compile the module:
%
% $ erl
%
% 1> c(fft).
% 2> fft:start().
%
% Or you can compile the program to Erlang bytecode and interpret it, using the following commands:
%
% $ erlc fft.erl
% $ erl -noshell -s fft start -s init stop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module attributes:
-module(anyfft).
-export([ start/0, recursive_fft/1, split_recursive_fft/1 ]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import needed modules:
-import(io, [ format/2, fwrite/1, fwrite/2 ]).                 % Input and output;
-import(math, [ atan2/2, cos/1, sin/1, pi/0, sqrt/1, pow/2 ]). % Math operations;
-import(lists, [ foldl/3, map/2, nth/2, seq/2, seq/3, zip/2, zipwith/3, duplicate/2, flatten/1 ]).     % List operations;
-import(timer, [ tc/2 ]).                                      % Time measurements;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mini-library to deal with complex numbers.
-type complex()::{ float(), float() }.

-spec csum(Z1::complex(), Z2::complex()) -> complex().         % Add two complex numbers;
csum({R1, I1}, {R2, I2}) -> { R1+R2, I1+I2 }.

-spec csub(Z1::complex(), Z2::complex()) -> complex().         % Subtract two complex numbers;
csub({R1, I1}, {R2, I2}) -> { R1-R2, I1-I2 }.

-spec cmul(Z1::complex(), Z2::complex()) -> complex().         % Multiply two complex numbers;
cmul({R1, I1}, {R2, I2}) -> { R1*R2 - I1*I2, R1*I2 + I1*R2 }.

-spec to_polar(Z::complex()) -> complex().                     % Convert to polar coordinates;
to_polar({R, I}) -> { sqrt(R*R + I*I), atan2(I, R) }.

-spec to_rect(Z::complex()) -> complex().                      % Convert to rectangular coordinates;
to_rect({Re, Th}) -> { Re*cos(Th), Re*sin(Th) }.

-spec cpow(Z::complex(), A::float()) -> complex().             % Power by a number;
cpow(Z, A) ->
    { Re, Th } = to_polar(Z),
    to_rect({ pow(Re, A), A*Th }).

-spec cexp(A :: float()) -> complex().                         % Complex exponential of an angle;
cexp(A) -> { cos(A), sin(A) }.

-spec c_to_string(Z :: complex()) -> string.
c_to_string(Z) -> format("(~8.4f, ~8.4f)", tuple_to_list(Z)).

% These versions of the functions return anonymous functions that will simplify the expressions in
% the functions below:
csum() -> fun(Z1, Z2) -> csum(Z1, Z2) end.             % Anonymous function to the sum;
csub() -> fun(Z1, Z2) -> csub(Z1, Z2) end.             % Anonymous function to the difference;
cmul() -> fun(Z1, Z2) -> cmul(Z1, Z2) end.             % Anonymous function to the product;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pretty printing of an array of complex numbers, used to inspect results.
%
% Parameters:
%   X
%     A vector of complex numbers, according to the definition above;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
complex_show([]) -> ok;
complex_show([H|T]) ->
    fwrite("~s~n", [ c_to_string(H) ]),
    complex_show(T).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure execution time through repeated calls to a (Fast) Fourier Transform function.
%
% Parameters:
%   F
%     Function to be called;
%   Size
%     Power of two of the size of the vector on which the transform will be applied;
%   Repeats
%     Number of times the function will be called. Defaults to REPEAT.
%
% Returns:
%   The average execution time for that function with a vector of the given size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec time_it(F::function(), Size::integer(), Repeats::integer()) -> float().
time_it(F, Size, Repeats) ->
    X = [ { J, 0.0 } || J <- seq(0, Size-1) ],         % Generate a vector;
    time_it(F, X, Repeats, 0.0) / (Repeats*1.0e6).     % Repeated calls and average;

time_it(_, _, 0, T) -> T;
time_it(F, X, Repeats, T) ->
    { T0, _ } = tc(F, [ X ]),                          % Time one iteration;
    time_it(F, X, Repeats-1, T+T0).                    % Sums with time from other iterations;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
% This implementation uses list functions to compute.
%
% This version is a tail-recursive implementation. It has the same parameters and behaviour of the
% non-tail recursive version, but is probably more time and memory efficient.
%
% Parameters:
%   x
%     The vector of which the DFT will be computed. Given the nature of the implementation, there is
%     no restriction on the size of the vector, although it will almost always be called with a
%     power of two size to give a fair comparison;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec direct_ft(X::[ complex() ]) -> [ complex() ].
direct_ft(X) -> direct_ft(X, 0, [ ]).

direct_ft(X, K, TX) when K == length(X) -> TX;
direct_ft(X, K, TX) ->
    Nx    = length(X),                                         % Length of the vector;
    W     = cexp(-2*pi()/Nx),                                  % Twiddle factor;
    Wk    = map(fun(N)->cpow(W, K*N) end, seq(0, Nx-1)),
    XnWkn = zipwith(cmul(), X, Wk),                            % Each term of the summation;
    TXk   = foldl(csum(), { 0.0, 0.0 }, XnWkn),                % The value of the kth component;
    direct_ft(X, K+1, TX ++ [ TXk ]).                          % Recursive call;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smallest prime factor of a given number. If the argument is prime itself, then it is the return
% value.
%
% Parameters:
%   N
%     Number to be inspected;
%   K
%     The prime number being inspected at the moment (control variable).
%
% Returns:
%   The smallest prime factor, or the number itself if it is already a prime.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec factor(N::integer(), K::integer()) -> integer().
factor(N, K) ->
    if
        N rem K == 0 -> K;                     % When remainder is zero, factor is found;
        K == N div 2 -> N;                     % Search up to half the number;
        true -> factor(N, K+1)
    end.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicate the elements of a list a given number of times.
%
% Parameters:
%  X
%    List to be replicated;
%  Repeats
%    Number of replications.
%
% Returns:
%   One list, with `Repeats` replications of the elements of the original list, in order.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec replicate(X::[ any() ], Repeats::integer()) -> [ any() ].
replicate(X, Repeats) -> flatten(duplicate(trunc(Repeats), X)).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
% complexity. This version uses comprehension lists to split vectors in even and odd samples.
%
% Parameters:
%   x
%     The vector of which the DFT will be computed. This should have length power of 2;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec recursive_fft(X::[ complex() ]) -> [ complex() ].
recursive_fft(X) ->
    Nx = length(X),
    N1 = factor(Nx, 2),
    if
        Nx == N1 -> direct_ft(X);              % Defer to the Direct FT;
        true -> recursive_fft(X, Nx, N1, 0)    % Summation loop;
    end.

recursive_fft(_, Nx, N1, J) when J == N1 -> duplicate(Nx, { 0.0, 0.0 });
recursive_fft(X, Nx, N1, J) ->
    Xj = recursive_fft([ nth(K, X) || K <- seq(J+1, Nx, N1) ]),
    Wkj = map(fun(K) -> cexp(-2*pi()*K*J/Nx) end, seq(0, Nx-1)),
    XjWkj = zipwith(cmul(), replicate(Xj, N1), Wkj),
    zipwith(csum(), recursive_fft(X, Nx, N1, J+1), XjWkj).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split the sequence X in N1 sequences of the same length. The length of X must be an integer
% multiple of N1.
%
% Parameters:
%   X
%     The vector to be processed;
%   N1
%     The number of sequences to be split;
%   J
%     Index of the sequence being split at the moment (control variable).
%
% Returns:
%   A list of complex-number vectors with the result of the splitting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec split(X::[ any() ], N1::integer(), J::integer()) -> [ [ any() ] ].
split(_, N1, J) when J > N1 -> [ ];
split(X, N1, J) ->
    Nx = length(X),
    Xj = [ nth(N, X) || N <- seq(J, Nx, N1) ],
    [ Xj | split(X, N1, J+1) ].


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
% complexity. This version uses an external function to split vectors in even and odd samples.
%
% Parameters:
%   x
%     The vector of which the DFT will be computed. This should have length power of 2;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec split_recursive_fft(X::[ complex() ]) -> [ complex() ].
split_recursive_fft(X) ->
    Nx = length(X),
    N1 = factor(Nx, 2),
    if
        Nx == N1 -> direct_ft(X);
        true ->
            Xjs = split(X, N1, 1),
            split_recursive_fft(Xjs, Nx, N1, 0)
    end.

split_recursive_fft([ ], Nx, _, _) -> duplicate(Nx, { 0.0, 0.0 });
split_recursive_fft([ Xj | Xjs ], Nx, N1, J) ->
    TXj = split_recursive_fft(Xj),
    Wkj = map(fun(K) -> cexp(-2*pi()*K*J/Nx) end, seq(0, Nx-1)),
    XjWkj = zipwith(cmul(), replicate(TXj, N1), Wkj),
    zipwith(csum(), split_recursive_fft(Xjs, Nx, N1, J+1), XjWkj).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_loop([], _) when R > RMax -> ok;
main_loop([N|T], Repeats) ->
    DTime  = time_it(fun(X) -> direct_ft(X) end, N, Repeats),
    RTime  = time_it(fun(X) -> recursive_fft(X) end, N, Repeats),
    STime  = time_it(fun(X) -> split_recursive_fft(X) end, N, Repeats),
    fwrite("| ~7w | ~7w | ~7.4f | ~7.4f | ~7.4f |~n", [ N, N*N, DTime, RTime, STime ]),
    main_loop(T, Repeats).


start() ->

    Sizees = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ]

    % Start by printing the table with time comparisons:
    fwrite("+---------+---------+---------+---------+---------+~n"),
    fwrite("|    N    |   N^2   | Direct  | Recurs. | SpRec.  |~n"),
    fwrite("+---------+---------+---------+---------+---------+~n"),

    % Try it with vectors with size ranging from 32 to 1024 samples:
    main_loop(Sizes, 50),

    fwrite("+---------+---------+---------+---------+---------+~n").
