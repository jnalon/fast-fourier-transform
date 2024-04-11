function Y = vec_recursive_nfft(X)

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
    N1 = prime_factor(N);                      % Find the smallest factor of the vector length;
    if N1 == N                                 % If the length is prime itself,
        Y = direct_ft(X);                      %    the transform is given by the direct form;
    else
        N2 = N / N1;                           % Decompose in two factors, N1 being prime;
        Y = zeros(1, N);                       % Accumulate the results;
        k = 0:N-1;
        Wk = exp(-2*i*pi*k/N);                 % Twiddle factors;
        Wkj = ones(1, N);
        for j = 1:N1                                   % Compute every subsequence of size N2;
            Xj = recursive_nfft(X(j:N1:end));          % Recursively compute the Fourier Transform;
            k2 = mod(k, N2) + 1;
            Y = Y + Xj(k2) .* Wkj;             % Recombine results;
            Wkj = Wkj .* Wk;                   % Update twiddle factors;
        end
    end

end
