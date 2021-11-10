classdef SymmetricGroup < replab.SignedPermutationGroup
% Describes the signed permutation group over {-n...-1, 1...n} where n = domainSize


    methods

        function self = SymmetricGroup(domainSize, generators, type, niceIsomorphism)
            order = replab.util.factorial(domainSize)*replab.util.multiplyIntegers(ones(1, domainSize)*2);
            self@replab.SignedPermutationGroup(domainSize, generators, 'niceIsomorphism', niceIsomorphism, 'type', type, 'order', order);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('Signed permutations acting on {-%d..-1 1..%d}', self.domainSize, self.domainSize);
        end

        % Domain

        function s = sample(self)
            n = self.domainSize;
            s = randperm(n) .* (randi([0 1], 1, n)*2-1);
        end

        % FiniteGroup

        function b = contains(self, g)
            assert(length(g) == self.domainSize, 'Signed permutation in wrong domain');
            assert(all(g ~= 0) && all(abs(g) <= self.domainSize), 'Signed permutation has out of range coefficients');
            b = true;
        end

    end

end
