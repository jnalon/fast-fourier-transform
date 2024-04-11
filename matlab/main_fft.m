%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- MATLAB/Octave Version
% This version implements Cooley-Tukey algorithm for powers of 2 only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was tested with Octave. All you need to do to run this program is to invoque the
% interpreter:
%
% $ octave main_fft.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions:
REPEAT = 50;                                   % Number of executions to compute average time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start by printing the table with time comparisons:
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');
fprintf('|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  | Intern. |\n');
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');

% Try it with vectors with size ranging from 32 to 1024 samples:
for r = 5:10

    % Compute the average execution time:
    n = 2^r;
    dtime = time_it(@direct_ft, n, REPEAT);
    rtime = time_it(@recursive_fft, n, REPEAT);
    itime = time_it(@iterative_fft, n, REPEAT);
    ptime = time_it(@fft, n, REPEAT);

    % Print the results:
    fprintf('| %7d | %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n', n, n*n, r*n, dtime, rtime, itime, ptime);

end
fprintf('+---------+---------+---------+---------+---------+---------+---------+\n');

