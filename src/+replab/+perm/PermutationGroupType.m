classdef PermutationGroupType < replab.FiniteGroupType
% The group type for permutation over a given domain size

    methods (Static)

        function G = make(domainSize)
        % Constructs the permutation group type for a given domain size
        %
        % This static method keeps the constructed copies in cache.
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
        %
        % Returns:
        %   `.PermutationGroupType`: The constructed or cached permutation group type
            persistent cache
            if isempty(cache)
                cache = cell(1, 0);
            end
            if domainSize+1 > length(cache) || isempty(cache{domainSize+1})
                cache{1,domainSize+1} = replab.perm.PermutationGroupType(domainSize);
            end
            G = cache{domainSize+1};
        end

    end

    properties (SetAccess = protected)
        domainSize % (integer): The integer ``d``, as this group acts on ``{1, ..., d}``
    end

    methods (Access = protected) % Constructor

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
            s = sprintf('Type: Permutations acting on 1..%d', self.domainSize);
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

        function c = compare(self, x, y)
            v = x - y;
            ind = find(v ~= 0, 1);
            c = [sign(v(ind)) 0];
            c = c(1);
        end

        function G = groupWithGenerators(self, generators, varargin)
            G = replab.PermutationGroup(self.domainSize, generators, varargin{:});
        end

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.perm.PermutationGroupType') && self.domainSize == otherType.domainSize;
        end

    end

end
