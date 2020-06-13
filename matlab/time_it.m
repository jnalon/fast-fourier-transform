function t = time_it(f, size, repeat)

% Auxiliary function: T = TIME_IT(F, SIZE, REPEAT)
% Measure execution time through repeated calls to a (Fast) Fourier Transform function.
%
% Parameters:
%   F
%     Function to be called;
%   SIZE
%     Number of elements of the vector on which the transform will be applied;
%   REPEAT
%     Number of times the function will be called. Defaults to REPEAT.
%
% Returns:
%   The average execution time for that function with a vector of the given size.Computa o tempo
%   médio de execução de N repetições de uma função.

    x = 0:size-1;                      % Generate a vector;
    t0 = cputime;                      % Start a timer;
    for j = 1:repeat,                  % Repeated calls;
        y = f(x);
    end
    t = (cputime - t0) / repeat;       % Compute average;

end
