classdef OrthogonalGroup < replab.CompactGroup
% Describes the group of n x n orthonormal (real) matrices

    properties
        n % integer: Dimension of the orthogonal group
        sparse % (logical): Whether to preserve sparse matrices
    end

    properties (Access = protected)
        parent % replab.domain.Matrices: Domain of square real matrices
    end

    methods

        function self = OrthogonalGroup(n, sparse)
            self.n = n;
            self.parent = replab.domain.Matrices('R', n, n);
            if n < 2
                sparse = false;
            end
            if sparse
                self.identity = speye(n);
            else
                self.identity = eye(n);
            end
        end

        %% Str methods

        function s = headerStr(self)
            s = sprintf('%d x %d orthonormal matrices', self.n, self.n);
        end

        %% Domain methods

        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end

        function X = sample(self)
            X = self.sampleUniformly;
        end

        %% Monoid methods

        function Z = compose(self, X, Y)
            Z = X * Y;
        end

        %% Group methods

        function XInv = inverse(self, X)
            XInv = X';
        end

        %% CompactGroup methods

        function X = sampleUniformly(self)
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            X = self.parent.sample;
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

        %% Representations

        function rep = definingRep(self)
        % Returns the natural representation of this orthonormal group
            ii = replab.irreducible.Info([], []);
            rep = replab.Rep.lambda(self, 'R', self.n, true, ii, @(o) o, @(o) o');
        end

    end

end
