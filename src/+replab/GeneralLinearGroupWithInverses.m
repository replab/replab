classdef GeneralLinearGroupWithInverses < replab.Group & replab.domain.VectorSpace
% Describes the group of n x n invertible real or complex matrices, elements store inverses

    properties
        n % (integer): Size of the square matrices
    end

    properties (Access = protected)
        parent_ % (+replab.+domain.Matrices): General, not necessarily invertible matrices
    end

    methods

        function self = GeneralLinearGroupWithInverses(field, n)
            self.field = field;
            self.n = n;
            self.parent_ = replab.domain.Matrices(field, n, n);
            self.identity = [eye(n) eye(n)];
        end

        % Str

        function s = headerStr(self)
            s = sprintf('%d x %d invertible matrices in %s', self.n, self.n, self.field);
        end

        % Domain

        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X(:,1:self.n), Y(:,1:self.n));
        end

        function X = sample(self)
            X = self.parent_.sample;
            X = [X inv(X)];
            % a generic gaussian matrix is almost always invertible
        end

        % Semigroup/monoid/group

        function Z = compose(self, X, Y)
            n = self.n;
            Xinv = X(:,n+1:2*n);
            X = X(:,1:n);
            Yinv = Y(:,n+1:2*n);
            Y = Y(:,1:n);
            Z = [X*Y Yinv*Xinv];
        end

        function Xinv = inverse(self, X)
            n = self.n;
            Xinv = X(:,n+1:2*n);
            X = X(:,1:n);
            Xinv = [Xinv X];
        end

    end

end
