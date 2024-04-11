function t = test_it(f, size)

% Auxiliary function: T = TEST_IT(F, SIZE)
% Pretty print of input and output of the Fourier transform for visual inspection.
%
% Parameters:
%   F
%     Function to be called;
%   SIZE
%     Number of elements of the vector on which the transform will be applied;

    x = 0:size-1;                      % Generate a vector;
    X = f(x);
    fprintf('N = %d | Input | Output\n', size);
    for i = 1:size,
        fprintf('  %2d | (%8.4f, %8.4f) | (%8.4f, %8.4f)\n',
                i, real(x(i)), imag(x(i)), real(X(i)), imag(X(i)));
    end
    fprintf('------------------------------\n');

end
