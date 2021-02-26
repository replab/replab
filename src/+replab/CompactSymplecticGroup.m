classdef CompactSymplecticGroup < replab.CompactGroup
% Describes the group Sp(n) of n x n unitary quaternion matrices

    properties (SetAccess = protected)
        n % (integer): Dimension of group
    end

    methods

        function self = CompactSymplecticGroup(n)
            self.n = n;
            self.identity = replab.H(eye(n));
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('Unitary symplectic group Sp(%d)', self.n);
        end

        % Domain

        function b = eqv(self, X, Y)
            b = full(all(all(X == Y)));
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
            X = replab.numerical.randomUnitaryOver(self.n, 'H');
        end

        function b = hasReconstruction(self)
            b = true;
        end

        function [mu, R] = reconstruction(self)
            R = replab.SetProduct.identity(self);
            T = replab.TorusGroup(self.n);
            mu = T.morphismByFunction(self, @(x) replab.H(T.toMatrix(x)));
        end

    end

    methods % Representations

        function rep = definingRep(self, field)
        % Returns the defining representation of this orthogonal group
        %
        % The defining representation of $Sp(d)$ corresponds to a ``2d x 2d`` complex unitary matrix
        % or a ``4d x 4d`` real orthogonal matrix.
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
