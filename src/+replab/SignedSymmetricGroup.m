classdef SignedSymmetricGroup < replab.SignedPermutationGroup
% Describes the signed permutation group over {-n...-1, 1...n} where n = domainSize

    methods (Static)

        function G = make(n)
        % Constructs the group of all signed permutation over a given domain size
        %
        % This static method keeps the constructed copies in cache.
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
        %
        % Returns:
        %   `.SignedSymmetricGroup`: The constructed or cached symmetric group
            persistent cache
            if isempty(cache)
                cache = cell(1, 0);
            end
            if n+1 > length(cache) || isempty(cache{n+1})
                cache{1,n+1} = replab.SignedSymmetricGroup(n);
            end
            G = cache{n+1};
        end

    end

    methods (Access = protected)

        function self = SignedSymmetricGroup(domainSize)
        % Constructs the group of all signed permutations over a given domain size
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
            generators = {[-1 2:domainSize]};
            if domainSize > 1
                generators{1,end+1} = [2:domainSize 1];
            end
            if domainSize > 2
                generators{1,end+1} = [2 1 3:domainSize];
            end
            type = 'self';
            ds = self.domainSize;
            o = replab.util.factorial(ds)*replab.util.multiplyIntegers(ones(1, ds)*2);
            % faster than o = factorial(vpi(domainSize))*vpi(2)^domainSize;
            self@replab.SignedPermutationGroup(domainSize, generators, 'type', type);
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
