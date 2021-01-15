classdef OrthogonalGroup < replab.CompactGroup
% Describes the group of n x n orthonormal (real) matrices

    properties
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
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            X = randn(self.n, self.n);
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

    end

    methods % Representations

        function rep = definingRep(self)
        % Returns the defining representation of this orthogonal group
        %
        % The defining representation of $O(d)$ corresponds to a ``d x d`` orthogonal matrix.
        %
        % Returns:
        %   `.Rep`: Defining representation
            rep = replab.rep.DefiningRep(self);
        end

    end

end
