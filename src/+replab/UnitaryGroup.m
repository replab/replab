classdef UnitaryGroup < replab.CompactGroup
% Describes the group of n x n unitary (complex) matrices

    properties (SetAccess = protected)
        n % (integer): Dimension of the unitary group
    end

    properties (Access = protected)
        parent % (`+replab.+domain.Matrices`): Domain of square complex matrices
    end

    methods

        function self = UnitaryGroup(n, identity)
        % Constructs the unitary group
        %
        % Args:
        %   n (integer): Size of the ``n x n`` matrices
        %   sparse (logical, optional): Whether to preserve sparse matrices, default value false
            self.n = n;
            self.parent = replab.domain.Matrices('C', n, n);
            if nargin < 2
                identity = eye(n);
            end
            self.identity = identity;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('%d x %d unitary matrices', self.n, self.n);
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
            X = replab.numerical.randomUnitaryOver(self.n, 'C');
        end

        function b = hasReconstruction(self)
            b = true;
        end

        function [R, mu] = reconstruction(self)
            R = replab.SetProduct.identity(self);
            T = replab.TorusGroup(self.n);
            mu = T.morphismByFunction(self, @(x) T.toMatrix(x));
        end

    end

    methods % Representations

        function rep = definingRep(self, field)
        % Returns the defining representation of this unitary group
        %
        % The defining representation of $U(d)$ corresponds to a ``d x d`` unitary matrix.
        %
        % Args:
        %   field ('R', 'C', optional): Field over which to define the matrices, default: C
        %
        % Returns:
        %   `.Rep`: Defining representation
            if nargin < 2 || isempty(field)
                field = 'C';
            end
            rep = replab.rep.DefiningRep(self, field);
        end

    end

end
