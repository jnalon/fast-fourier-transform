%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier Transform -- MATLAB/Octave Version
% This program performs a simple test of the implementations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version was tested with Octave. All you need to do to run this program is to invoque the
% interpreter:
%
% $ octave main_test.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions:
REPEAT = 50;                                   % Number of executions to compute average time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for the implementations of the DFT of power-of-2 length:
fprintf('Direct FT - \n');
test_it(@direct_ft, 8);
fprintf('Recursive FT - \n');
test_it(@recursive_fft, 8);
fprintf('Iterative FT - \n');
test_it(@iterative_fft, 8);
fprintf('Direct FT - \n');
test_it(@direct_ft, 16);
fprintf('Recursive FT - \n');
test_it(@recursive_fft, 16);
fprintf('Iterative FT - \n');
test_it(@iterative_fft, 16);

% Tests for the implementations of the DFT of composite length:
fprintf('Direct FT - \n');
test_it(@direct_ft, 12);
fprintf('Recursive FT - \n');
test_it(@recursive_nfft, 12);
fprintf('Vector Recursive FT - \n');
test_it(@vec_recursive_nfft, 12);
