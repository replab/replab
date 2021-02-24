classdef OrthogonalGroup < replab.CompactGroup
% Describes the group of n x n orthonormal (real) matrices

    properties (SetAccess = protected)
        n % (integer): Dimension of the orthogonal group
    end

    properties (Access = protected)
        parent % (`+replab.+domain.Matrices`): Domain of square real matrices
    end

    methods

        function self = OrthogonalGroup(n, identity)
            self.n = n;
            self.parent = replab.domain.Matrices('R', n, n);
            if nargin < 2
                identity = eye(n);
            end
            self.identity = identity;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('%d x %d orthonormal matrices', self.n, self.n);
        end

        % Domain

        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end

        % Monoid

        function Z = compose(self, X, Y)
            Z = X * Y;
        end

        % Group

        function XInv = inverse(self, X)
            XInv = X';
        end

        % CompactGroup

        function X = sample(self)
            X = replab.numerical.randomUnitaryOver(self.n, 'R');
        end

    end

    methods % Representations

        function rep = definingRep(self, field)
        % Returns the defining representation of this orthogonal group
        %
        % The defining representation of $O(d)$ corresponds to a ``d x d`` orthogonal matrix.
        %
        % Args:
        %   field ('R', 'C'): Field over which to define the representation
        %
        % Returns:
        %   `.Rep`: Defining representation
            if nargin < 2 || isempty(field)
                field = 'R';
            end
            rep = replab.rep.DefiningRep(self, field);
        end

    end

end
