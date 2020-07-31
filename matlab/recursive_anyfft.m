function Y = recursive_anyfft(X)

% Y = RECURSIVE_ANYFFT(X)
% Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller complexity
% than the direct FT, though the exact value is difficult to compute.
%
% Parameters:
%   X
%     The vector of which the FFT will be computed. Its length must be a composite number, or else
%     the computation will be defered to the direct FT, and there will be no efficiency gain.
%
% Returns:
%    A complex-number vector of the same size, with the coefficients of the DFT.

    N = length(X);                             % Length of the vector;
    N1 = factor(N);                            % Find the smallest factor of the vector length;
    if N1 == N                                 % If the length is prime itself,
        Y = direct_ft(X);                      %    the transform is given by the direct form;
    else
        N2 = N / N1;                           % Decompose in two factors, N1 being prime;
        Y = zeros(1, N);                       % Accumulate the results;
        W = exp(-2*i*pi/N);                    % Twiddle factors;
        Wj = 1;
        for j = 1:N1                           % Compute every subsequence of size N2;
            Xj = recursive_anyfft(X(j:N1:end));
            Wkj = 1;
            for k = 1:N
                k2 = mod(k-1, N2) + 1;
                Y(k) = Y(k) + Xj(k2) * Wkj;    % Recombine results;
                Wkj = Wkj * Wj;                % Update twiddle factors;
            end
            Wj = Wj * W;
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = factor(n)

% K = FACTOR(N)
% Smallest prime factor of a given number. If the argument is prime itself, then it is the return
% value.
%
% Parameters:
%   N
%     Number to be inspected.
%
% Returns:
%  The smallest prime factor, or the number itself if it is already a prime.

    rn = floor(sqrt(n));                       % Search up to the square root of the number;
    for i = 2:rn
        if mod(n, i) == 0                      % When remainder is zero, factor is found;
            k = i;
            return
        end
    end
    k = n;

end
