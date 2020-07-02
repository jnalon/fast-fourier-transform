%% -*- coding: utf-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- Python 3 Version
% This version implements Cooley-Tukey algorithm for powers of 2 only.
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
-module(fft).
-export([ start/0 ]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import needed modules:
-import(io, [ format/2, fwrite/1, fwrite/2 ]).                 % Input and output;
-import(math, [ atan2/2, cos/1, sin/1, pi/0, sqrt/1, pow/2 ]). % Math operations;
-import(lists, [ foldl/3, map/2, nth/2, seq/2, seq/3, zip/2, zipwith/3 ]).     % List operations;
-import(timer, [ tc/2 ]).                                      % Time measurements;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mini-library to deal with complex numbers.
-type complex() :: { float(), float() }.

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

-spec cpow(Z::complex(), A :: float()) -> complex().           % Power by a number;
cpow(Z, A) ->
    { Re, Th } = to_polar(Z),
    to_rect({ pow(Re, A), A*Th }).

-spec cexp(A::float()) -> complex().                           % Complex exponential of an angle;
cexp(A) -> { cos(A), sin(A) }.

-spec c_to_string(Z :: complex()) -> string.
c_to_string(Z) -> format("(~8.4f, ~8.4f)", tuple_to_list(Z)).

% These versions of the functions return anonymous functions that will simplify the expressions in
% the functions below:
csum() -> fun(Z1, Z2) -> csum(Z1, Z2) end.                     % Anonymous function to the sum;
csub() -> fun(Z1, Z2) -> csub(Z1, Z2) end.                     % Anonymous function to the difference;
cmul() -> fun(Z1, Z2) -> cmul(Z1, Z2) end.                     % Anonymous function to the product;


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
% This implementation uses list comprehensions to compute.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. Given the nature of the implementation, there is
%     no restriction on the size of the vector, although it will almost always be called with a
%     power of two size to give a fair comparison;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec direct_ft(X::[ complex() ]) -> [ complex() ].
direct_ft(X) -> direct_ft(X, 0).

direct_ft(X, K) when K == length(X) -> [ ];                    % FT of an empty vector;
direct_ft(X, K) ->                                             % Computes kth component;
    Nx    = length(X),                                         % Length of the vector;
    W     = cexp(-2*pi()/Nx),                                  % Twiddle factor;
    Wk    = [ cpow(W, K*N)  || N <- seq(0, Nx-1) ],
    XnWkn = [ cmul(Xn, Wkn) || { Xn, Wkn } <- zip(X, Wk) ],    % Each term of the summation;
    TXk   = foldl(csum(), { 0.0, 0.0 }, XnWkn),                % The value of the kth component;
    [ TXk | direct_ft(X, K+1) ].                               % Recursive call;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
% This implementation uses list functions to compute.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. Given the nature of the implementation, there is
%     no restriction on the size of the vector, although it will almost always be called with a
%     power of two size to give a fair comparison;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec list_direct_ft(X::[ complex() ]) -> [ complex() ].
list_direct_ft(X) -> list_direct_ft(X, 0).

list_direct_ft(X, K) when K == length(X) -> [ ];               % FT of an empty vector;
list_direct_ft(X, K) ->                                        % Computes kth component;
    Nx    = length(X),                                         % Length of the vector;
    W     = cexp(-2*pi()/Nx),                                  % Twiddle factor;
    Wk    = map(fun(N)->cpow(W, K*N) end, seq(0, Nx-1)),
    XnWkn = zipwith(cmul(), X, Wk),                            % Each term of the summation;
    TXk   = foldl(csum(), { 0.0, 0.0 }, XnWkn),                % The value of the kth component;
    [ TXk | list_direct_ft(X, K+1) ].                          % Recursive call;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
% This implementation uses list comprehensions to compute.
%
% This version is a tail-recursive implementation. It has the same parameters and behaviour of the
% non-tail recursive version, but is probably more time and memory efficient.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. Given the nature of the implementation, there is
%     no restriction on the size of the vector, although it will almost always be called with a
%     power of two size to give a fair comparison;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec tr_direct_ft(X::[ complex() ]) -> [ complex() ].
tr_direct_ft(X) -> tr_direct_ft(X, 0, [ ]).

tr_direct_ft(X, K, TX) when K == length(X) -> TX;
tr_direct_ft(X, K, TX) ->
    Nx    = length(X),                                         % Length of the vector;
    W     = cexp(-2*pi()/Nx),                                  % Twiddle factor;
    Wk    = [ cpow(W, K*N)  || N <- seq(0, Nx-1) ],
    XnWkn = [ cmul(Xn, Wkn) || { Xn, Wkn } <- zip(X, Wk) ],    % Each term of the summation;
    TXk   = foldl(csum(), { 0.0, 0.0 }, XnWkn),                % The value of the kth component;
    tr_direct_ft(X, K+1, TX ++ [ TXk ]).                       % Recursive call;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
% This implementation uses list functions to compute.
%
% This version is a tail-recursive implementation. It has the same parameters and behaviour of the
% non-tail recursive version, but is probably more time and memory efficient.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. Given the nature of the implementation, there is
%     no restriction on the size of the vector, although it will almost always be called with a
%     power of two size to give a fair comparison;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec list_tr_direct_ft(X::[ complex() ]) -> [ complex() ].
list_tr_direct_ft(X) -> list_tr_direct_ft(X, 0, [ ]).

list_tr_direct_ft(X, K, TX) when K == length(X) -> TX;
list_tr_direct_ft(X, K, TX) ->
    Nx    = length(X),                                         % Length of the vector;
    W     = cexp(-2*pi()/Nx),                                  % Twiddle factor;
    Wk    = map(fun(N)->cpow(W, K*N) end, seq(0, Nx-1)),
    XnWkn = zipwith(cmul(), X, Wk),                            % Each term of the summation;
    TXk   = foldl(csum(), { 0.0, 0.0 }, XnWkn),                % The value of the kth component;
    list_tr_direct_ft(X, K+1, TX ++ [ TXk ]).                  % Recursive call;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implementation uses list comprehensions, which are a more natural way to implement repeated
% computations in functional languages. Here is an explanation of every term:
% [ foldl(fun(Z1, Z2) -> csum(Z1, Z2) end,             % Summation;
%         { 0, 0 },                                    % Accumulator;
%         [
%           cmul(Xn,                                   %   of the nth sample
%                cexp(-2*pi()*K*N/length(X)))          %   times the twiddle factor
%           || { Xn, N } <- enumerate(X)               %   over the interval of samples;
%         ]
%         || K <- seq(0, length(X)-1) ]                % Repeat for every coefficient;
%
% Here is the analysis equation (in LaTeX form) for comparison:
%
%   X[k] = \sum_{0}^{N-1} x[n] e^{-j 2 \pi k n / N} \; k = 0 \ldots N-1
%
% This implementation is not, of course, very efficient, since it doesn't take advantage of the
% regularity of the twiddle factors, and computes complex exponentials for every term.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec enumerate(X::[ any() ]) -> [ { integer(), any() } ].
enumerate(X) -> zip(seq(0, length(X)-1), X).

-spec lc_dft(X :: [ complex() ]) -> [ complex() ].
lc_dft(X) -> [ foldl(csum(), { 0.0, 0.0 },
                     [ cmul(Xn, cexp(-2*pi()*K*N/length(X))) || { N, Xn } <- enumerate(X) ])
               || K <- seq(0, length(X)-1) ].


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
% complexity. This version uses comprehension lists to split vectors in even and odd samples.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. This should have length power of 2;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec recursive_fft(X::[ complex() ]) -> [ complex() ].
recursive_fft(X) when length(X) == 1 -> X;
recursive_fft(X) ->
    Nx  = length(X),
    Wk  = [ cexp(-2*pi()*K/Nx) || K <- seq(0, Nx div 2 - 1) ],
    Xe = [ nth(N, X) || N <- seq(1, Nx, 2) ],
    Xo = [ nth(N, X) || N <- seq(2, Nx, 2) ],
    TXe = recursive_fft(Xe),
    TXo = zipwith(cmul(), Wk, recursive_fft(Xo)),
    zipwith(csum(), TXe, TXo) ++ zipwith(csub(), TXe, TXo).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
% complexity. This version uses a external function to split vectors in even and odd samples.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. This should have length power of 2;
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec split(X::[ any() ]) -> { [ any() ], [ any() ] }.
split([]) -> { [], [] };
split(X) ->
    [ Xe | PTail ] = X,
    [ Xo | Tail ] = PTail,
    { XeTail, XoTail } = split(Tail),
    { [ Xe | XeTail ], [ Xo | XoTail ] }.

-spec split_recursive_fft(X::[ complex() ]) -> [ complex() ].
split_recursive_fft(X) when length(X) == 1 -> X;
split_recursive_fft(X) ->
    Nx  = length(X),
    Wk  = [ cexp(-2*pi()*K/Nx) || K <- seq(0, Nx div 2 - 1) ],
    { Xe, Xo } = split(X),
    TXe = split_recursive_fft(Xe),
    TXo = zipwith(cmul(), Wk, split_recursive_fft(Xo)),
    zipwith(csum(), TXe, TXo) ++ zipwith(csub(), TXe, TXo).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_loop(R, RMax, _) when R > RMax -> ok;
main_loop(R, RMax, Repeats) ->
    N = trunc(pow(2, R)),
    DTime  = time_it(fun(X) -> direct_ft(X) end, N, Repeats),
    LDTime = time_it(fun(X) -> list_direct_ft(X) end, N, Repeats),
    TTime  = time_it(fun(X) -> tr_direct_ft(X) end, N, Repeats),
    LTime  = time_it(fun(X) -> list_tr_direct_ft(X) end, N, Repeats),
    LCTime = time_it(fun(X) -> lc_dft(X) end, N, Repeats),
    RTime  = time_it(fun(X) -> recursive_fft(X) end, N, Repeats),
    STime  = time_it(fun(X) -> split_recursive_fft(X) end, N, Repeats),
    fwrite("| ~7w | ~7w | ~7w | ~7.4f | ~7.4f | ~7.4f | ~7.4f | ~7.4f | ~7.4f | ~7.4f |~n",
            [ N, N*N, N*R, DTime, LDTime, TTime, LTime, LCTime, RTime, STime ]),
    main_loop(R+1, RMax, Repeats).


start() ->

    % Start by printing the table with time comparisons:
    fwrite("+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+~n"),
    fwrite("|    N    |   N^2   | N logN  | Direct  | ListDir.| TRDir.  | LTRDir. | LCDir.  | Recurs. | Split R.|~n"),
    fwrite("+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+~n"),

    % Try it with vectors with size ranging from 32 to 1024 samples:
    main_loop(5, 10, 50),

    fwrite("+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+~n").
