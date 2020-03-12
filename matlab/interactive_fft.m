function Y = interactive_fft(X)

% Y = INTERACTIVE_FFT(X)
% Computes the Fast Fourier Ttransform using an interactive in-place decimation in time algorithm.
% This has O(N log_2(N)) complexity, and since there are less function calls, it will probably be
% marginally faster than the recursive versions.
%
% Parameters:
%   X
%     The vector of which the FFT will be computed. This should always be called with a vector of
%     a power of two length, or it will fail. No checks on this are made.
%
% Returns:
%  A complex-number vector of the same size, with the coefficients of the DFT.

    N = length(X);                             % Length of vector;
    r = floor(log(N)/log(2));                  % Number of bits;
    Y = X;                                     % Accumulates the results;
    for k = 0:N-1                              % Reorder the vector according to the
        l = bit_reverse(k, r);                 %   bit-reversed order;
        Y(l+1) = X(k+1);
    end
    step = 1;                                  % Auxililary for computation of twiddle factors;
    for k = 1:r,
        for l = 0:2*step:N-1
            W = exp(-i*pi/step);               % Twiddle factors;
            Wkn = 1;
            for n = 0:step-1,
                p = l + n + 1;
                q = p + step;
                Y(q) = Y(p) - Wkn*Y(q);        % Recombine results;
                Y(p) = 2*Y(p) - Y(q);
                Wkn = Wkn * W;                 % Update twiddle factors;
            end
        end
      step = 2*step;
   end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = bit_reverse(K, R)

% M = BIT_REVERSE(R, K)
% Computes the bit-reversed function of an integer number.
%
% Parameters:
%   K
%     The number to be bit-reversed;
%   R
%     The number of bits to take into consideration when reversing.
%
% Returns:
%   The number K, bit-reversed according to integers with R bits.

    L = 0;                                     % Accumulates the results;
    for i = 1:R,                               % Loop on every bit;
        L = 2*L;
        if mod(K, 2) == 1,                     % Tests less signficant bit and add;
            L = L + 1;
        end
        K = floor(K/2);                        % Tests next bit;
    end

end