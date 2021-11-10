classdef FiniteGroupType < replab.gen.FiniteGroupType
% A type for signed permutation groups

    methods (Static)

        function G = make(domainSize)
        % Constructs the signed permutation group type for a given domain size
        %
        % This static method keeps the constructed copies in cache
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
        %
        % Returns:
        %   `.FiniteGroupType`: The constructed or cached signed permutation group type
            persistent cache
            if isempty(cache)
                cache = cell(1, 0);
            end
            if domainSize+1 > length(cache) || isempty(cache{domainSize+1})
                cache{1,domainSize+1} = replab.signed.FiniteGroupType(domainSize);
            end
            G = cache{domainSize+1};
        end

    end

    properties (SetAccess = protected)
        domainSize % (integer): Domain size $d$ where elements act on $\{-d .. -1, 1 .. d\}$
        isomorphism % (`+replab.signed.NiceIsomorphism`): Standard isomorphism
    end

    methods (Access = protected)

        function self = FiniteGroupType(domainSize)
        % Constructs a signed permutation group type
        %
        % Args:
        %   domainSize (integer): Size of the domain
            self.identity = 1:domainSize;
            self.domainSize = domainSize;
            self.isomorphism = replab.signed.NiceIsomorphism(domainSize, self);
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(x == y);
        end

        function s = sample(self)
            n = self.domainSize;
            s = randperm(n) .* (randi([0 1], 1, n)*2-1);
        end

        % TotalOrder

        function c = compare(self, x, y)
            v = replab.SignedPermutation.toPermutation(x) - replab.SignedPermutation.toPermutation(y);
            ind = find(v ~= 0, 1);
            c = [sign(v(ind)) 0];
            c = c(1);
        end

        % Monoid

        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end

        % Group

        function y = inverse(self, x)
            n = self.domainSize;
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end

    end

end
