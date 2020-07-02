function Y = vec_recursive_anyfft(X)

% Y = VEC_RECURSIVE_ANYFFT(X)
% Fast Fourier Transform using a recursive decimation in time algorithm. This has smaller complexity
% than the direct FT, though the exact value is difficult to compute. In this implementation, loops
% are avoided by vectorizing the computation of the twiddle factors.
%
% Parameters:
%   X
%     The vector of which the FFT will be computed. It must be a composite number, or else the
%     computation will be defered to the direct FT, and there will be no efficiency gain.
%
% Returns:
%    A complex-number vector of the same size, with the coefficients of the DFT.

    N = length(X);                             % Length of the vector;
    N1 = factor(N);                            % Find the smallest factor of the vector length;
    if N1 == N                                 % If the length is prime itself,
        Y = direct_ft(X);                      %    the transform is given by the direct form;
    else
        N2 = N / N1;                                   % Decompose in two factors, N1 being prime;
        Y = zeros(1, N);                               % Accumulate the results;
        k = 0:N-1;
        for j = 1:N1                                   % Compute every subsequence of size N2;
            Xj = recursive_anyfft(X(j:N1:end));        % Recursively compute the Fourier Transform;
            Wkj = exp(-2*i*pi*k*j/N);
            k2 = mod(k, N2);
            Y = Y + Xj(k2+1) .* Wkj;           % Recombine results;
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
