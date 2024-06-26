%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- MATLAB/Octave Version
% This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was tested with Octave. All you need to do to run this program is to invoque the
% interpreter:
%
% $ octave main_anyfft.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions:
REPEAT = 50;                                   % Number of executions to compute average time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start by printing the table with time comparisons:
fprintf('+---------+---------+---------+---------+---------+---------+\n');
fprintf('|    N    |   N^2   | Direct  | Recurs. | VecRec. | Intern. |\n');
fprintf('+---------+---------+---------+---------+---------+---------+\n');

% Try it with vectors with the given sizes:
sizes = [ 2*3, 2*2*3, 2*3*3, 2*3*5, 2*2*3*3, 2*2*5*5, 2*3*5*7, 2*2*3*3*5*5 ];
for i = 1:length(sizes)

    % Compute the average execution time:
    n = sizes(i);
    dtime = time_it(@direct_ft, n, REPEAT);
    rtime = time_it(@recursive_nfft, n, REPEAT);
    vtime = time_it(@vec_recursive_nfft, n, REPEAT);
    ptime = time_it(@fft, n, REPEAT);

    % Print the results:
    fprintf('| %7d | %7d | %7.4f | %7.4f | %7.4f | %7.4f |\n', n, n*n, dtime, rtime, vtime, ptime);

end
fprintf('+---------+---------+---------+---------+---------+---------+\n');
