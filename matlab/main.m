%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- MATLAB/Octave Version
% This version implements Cooley-Tukey algorithm for powers of 2 only.
%
% José Alexandre Nalon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was tested with Octave. All you need to do to run this program is to invoque the
% interpreter:
%
% $ octave main.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitios:
REPEAT = 50;                                   % Number of executions to compute average time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starts by printing the table with time comparisons:
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');
fprintf('|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  | Interna |\n');
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');

% Try it with vectors with size ranging from 32 to 1024 samples:
for r = 5:10

    % Computes the average execution time:
    dtime = time_it(@direct_ft, r, REPEAT);
    rtime = time_it(@recursive_fft, r, REPEAT);
    itime = time_it(@interactive_fft, r, REPEAT);
    ptime = time_it(@fft, r, REPEAT);

    % Print the results:
    n = 2^r;
    fprintf('| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n', n, n*n, r*n, dtime, rtime, itime, ptime);

end
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');

