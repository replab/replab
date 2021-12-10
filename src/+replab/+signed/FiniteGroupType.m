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
            self.identity = 1:domainSize;
            self.domainSize = domainSize;
            n = self.domainSize;
            sourceGenerators = {[-1 2:n]};
            if n > 1
                sourceGenerators{1,end+1} = [2:n 1];
            end
            if n > 2
                sourceGenerators{1,end+1} = [2 1 3:n];
            end
            orderFun = @() replab.util.factorial(domainSize)*replab.util.multiplyIntegers(ones(1, domainSize)*2);
            sourceArgs = {'order', orderFun};
            targetType = replab.perm.PermutationGroupType.make(2*domainSize);
            self.finishConstruction(sourceGenerators, sourceArgs, targetType);
        end

    end

    methods % Implementations

        function S = makeParentGroup(self, generators, nice, niceIsomorphism)
            S = replab.signed.SymmetricGroup(self, generators, nice, niceIsomorphism);
        end

    end

    methods (Access = protected) % Implementations

        function G = groupFromNiceImage_(self, generators, nice, niceIsomorphism)
            G = replab.SignedPermutationGroup(self.domainSize, generators, 'type', self, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
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

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.signed.FiniteGroupType') && self.domainSize == otherType.domainSize;
        end

        % StaticFiniteGroupType

        function t = imageElement(self, s)
            t = replab.SignedPermutation.toPermutation(s);
        end

        function s = preimageElement(self, t)
            s = replab.SignedPermutation.fromPermutation(t);
        end

    end

end
