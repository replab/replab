classdef PermutationGroupType < replab.FiniteGroupType
% The group type for permutation over a given domain size

    properties (SetAccess = protected)
        domainSize % integer: The integer ``d``, as this group acts on ``{1, ..., d}``
    end

    methods % Constructor

        function self = PermutationGroupType(domainSize)
        % Constructs a permutation group type
        %
        % Args:
        %   domainSize (integer): Size of the domain
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = 'Finite group type of permutations';
        end

        % Domain

        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        function s = sample(self)
            s = randperm(self.domainSize);
        end

        % Monoid

        function z = compose(self, x, y)
            z = x(y);
        end

        % Group

        function y = inverse(self, x)
            n = length(x);
            y = zeros(1, n);
            y(x) = 1:n;
        end

        function z = composeWithInverse(self, x, y)
            z = zeros(1, length(x));
            z(y) = x;
        end

        function z = leftConjugate(self, x, y)
            z = zeros(1, length(x));
            % x y xInv
            z(x) = x(y);
        end

        % FiniteGroupType

        function G = groupWithGenerators(self, generators, varargin)
            G = replab.PermutationGroup(self, generators, varargin{:});
        end

    end

end
