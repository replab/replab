classdef FiniteGroupType < replab.gen.StaticFiniteGroupType
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
    end

    methods (Access = protected)

        function self = FiniteGroupType(domainSize)
        % Constructs a signed permutation group type
        %
        % Args:
        %   domainSize (integer): Size of the domain
            self.niceType = replab.PermutationGroupType.make(2*domainSize);
            self.identity = 1:domainSize;
            self.domainSize = domainSize;
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

        % FiniteGroupType

        function G = groupWithGenerators(self, generators, varargin)
            G = replab.SignedPermutationGroup(self.domainSize, generators, 'type', self, varargin{:});
        end

        % StaticFiniteGroupType

        function t = imageElement(self, s)
            t = replab.SignedPermutation.toPermutation(s);
        end

        function S = makeSource(self, generators, niceIsomorphism)
            S = replab.signed.SymmetricGroup(self.domainSize, generators, self, niceIsomorphism);
        end

        function s = preimageElement(self, t)
            s = replab.SignedPermutation.fromPermutation(t);
        end

        function G = sourceGenerators(self)
            n = self.domainSize;
            G = {[-1 2:n]};
            if n > 1
                G{1,end+1} = [2:n 1];
            end
            if n > 2
                G{1,end+1} = [2 1 3:n];
            end
        end

    end

end
