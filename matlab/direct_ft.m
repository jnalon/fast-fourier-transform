function Y = direct_ft(X)

% Y = DIRECT_FT(X)
% Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2) complexity.
%
% Parameters:
%   X
%     The vector of which the DFT will be computed. Given the nature of the implementation, there
%     is no restriction on the size of the vector, although it will almost always be called with
%     a power of two size to give a fair comparison.
%
% Returns:
%   A complex-number vector of the same size, with the coefficients of the DFT.

    N = length(X);                             % Length of the vector;
    Y = zeros(1, N);                           % Accumulate the results;
    W = exp(-2*i*pi/N);                        % Twiddle factors;
    Wk = 1;
    for k = 1:N                                % Compute the kth coefficient;
        Wkn = 1;
        for n = 1:N                            %   Operate the summation;
            Y(k) = Y(k) + X(n) .* Wkn;         %     Compute every term;
            Wkn = Wkn * Wk;                    % Update twiddle factors;
        end
        Wk = Wk * W;
    end

end
