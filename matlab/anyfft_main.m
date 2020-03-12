%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- MATLAB/Octave Version
% This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
%
% José Alexandre Nalon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was tested with Octave. All you need to do to run this program is to invoque the
% interpreter:
%
% $ octave anyfft_main.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitios:
REPEAT = 50;                                   % Number of executions to compute average time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starts by printing the table with time comparisons:
fprintf('+---------+---------+---------+---------+---------+\n');
fprintf('|    N    |   N^2   | Direct  | Recurs. | Interna |\n');
fprintf('+---------+---------+---------+---------+---------+\n');

% Try it with vectors with size ranging from 32 to 1024 samples:
sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ];
for i = 1:length(sizes)

    % Computes the average execution time:
    n = sizes(i);
    dtime = time_it(@direct_ft, n, REPEAT);
    rtime = time_it(@recursive_anyfft, n, REPEAT);
    ptime = time_it(@fft, n, REPEAT);

    % Print the results:
    fprintf('| %7d | %7d | %7.4f | %7.4f | %7.4f |\n', n, n*n, dtime, rtime, ptime);

end
fprintf('+---------+---------+---------+---------+---------+\n');
