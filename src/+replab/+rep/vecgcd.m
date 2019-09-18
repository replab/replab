function g = vecgcd(vec)
% Returns the greatest common divisor of a vector of integers
%
% Args:
%   vec (vector of integer): Integer to compute the GCD of
%
% Returns:
%   integer: The computed GCD
    assert(length(vec) > 0);
    g = vec(1);
    for i = 2:length(vec)
        g = gcd(g, vec(i));
    end
end
