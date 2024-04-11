function k = prime_factor(n)

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
