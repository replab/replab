function generators = generatorsPrimeResidues(n)
% Generators of the Galois group
%
% Follows the algorithm in GeneratorsPrimeResidues / GAP 4
%
% Example:
%   >>> replab.numerical.integer.generatorsPrimeResidues(1)
%       {}
%   >>> replab.numerical.integer.generatorsPrimeResidues(4*3)
%       {[7], [5]}
%   >>> replab.numerical.integer.generatorsPrimeResidues(8*9*5)
%       {[271, 181], [281], [217]}

% Args:
%   n (integer): Positive integer
%
% Returns:
%   cell(1,\*) of integer(1,\*): A list describing the generators of the group of prime residues
    if n == 1
        generators = {};
        return
    end

    [primes, exponents] = replab.numerical.integer.collected(factor(n));

    generators = {};
    for i = 1:length(primes)
        ppart = primes(i)^exponents(i);
        rest = n / ppart;
        [~, X] = replab.numerical.integer.xgcd([ppart rest]);
        gcd = X(2)*rest;
        if primes(i) == 2
            if mod(ppart, 8) == 0
                generators{i} = [mod(-2*gcd+1, n) mod(4*gcd+1, n)];
            else
                generators{i} = mod(-2*gcd+1, n);
            end
        else
            generators{i} = mod((replab.numerical.integer.generator(ppart) - 1) * gcd + 1, n);
        end
    end
end
