classdef SignedPermutationGroupType < replab.GenericFiniteGroupType
% A type for signed permutation groups

    properties (SetAccess = protected)
        domainSize % (integer): Domain size $d$ where elements act on $\{-d .. -1, 1 .. d\}$
        genericIsomorphism % (`.GenericIsomorphism`): Isomorphism for all groups of this type
    end

    methods

        function self = SignedPermutationGroupType(domainSize)
        % Constructs a signed permutation group type
        %
        % Args:
        %   domainSize (integer): Size of the domain
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

        % GenericFiniteGroupType

        function iso = genericIsomorphism(self, varargin)

        end

    end

end
