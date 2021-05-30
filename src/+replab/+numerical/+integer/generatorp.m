function g = generatorp(p)
% Returns a generator for the cyclic group of multiplication modulo p (p prime)
%
% Also called a primitive root of p.
%
% Args:
%   p (integer): A prime number
%
% Returns:
%   integer: The generator
%
% Taken from https://people.cs.kuleuven.be/~dirk.nuyens/fast-cbc/
% (C) <dirk.nuyens@cs.kuleuven.ac.be>
    if ~isprime(p)
        error('p = %d must be prime', p)
    end
    primef = unique(factor(p-1));
    g = 2;
    i = 1;
    while i <= length(primef)
        if replab.numerical.integer.powmod(g, (p-1)/primef(i), p) == 1
            g = g + 1;
            i = 0;
        end
        i = i + 1;
    end
end
