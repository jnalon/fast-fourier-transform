function Y = direct_ft(X)

% Y = DIRECT_FT(X)
% Computes the Discrete Fourier Ttransform directly from the definition, an algorithm that has
% O(N^2) complexity.
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
    Y = zeros(1, N);                           % Accumulates the results;
    W = exp(-2*i*pi/N);                        % Twiddle factors;
    Wk = 1;
    for k = 1:N                                % Computes the kth coefficient;
        Wkn = 1;
        for n = 1:N                            %   Operates the summation;
            Y(k) = Y(k) + X(n) .* Wkn;         %     Computes every term;
            Wkn = Wkn * Wk;                    % Update twiddle factors;
        end
        Wk = Wk * W;
    end

end
