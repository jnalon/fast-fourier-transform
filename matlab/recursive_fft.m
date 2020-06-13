function Y = recursive_fft(X)

% Y = RECURSIVE_FFT(X)
% Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
% complexity.
%
% Parameters:
%   X
%     The vector of which the FFT will be computed. This should always be called with a vector of
%     a power of two length, or it will fail. No checks on this are made.
%
% Returns:
%    A complex-number vector of the same size, with the coefficients of the DFT.

    N = length(X);                             % Length of the vector;
    if N == 1                                  % A length-1 vector is its own FT;
        Y = X;
    else
        Y = zeros(1, N);                       % Length of the vector;
        Xe = recursive_fft(X(1:2:end));        % Transform of even samples;
        Xo = recursive_fft(X(2:2:end));        % Transform of odd samples;
        W = exp(-2*i*pi*(0:N/2-1)/N);          % Twiddle factors;
        WXo = W .* Xo;                         % Repeated computation;
        Y = [ Xe + WXo, Xe - WXo ];            % Recombine results;
    end

end
