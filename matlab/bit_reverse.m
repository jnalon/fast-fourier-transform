function L = bit_reverse(K, R)

% M = BIT_REVERSE(K, R)
% Bit-reversed version of an integer number.
%
% Parameters:
%   K
%     The number to be bit-reversed;
%   R
%     The number of bits to take into consideration when reversing.
%
% Returns:
%   The number K, bit-reversed according to integers with R bits.

    L = 0;                                     % Accumulate the results;
    for i = 1:R,                               % Loop on every bit;
        L = 2*L + mod(K, 2);
        K = floor(K/2);                        % Test next bit;
    end

end
