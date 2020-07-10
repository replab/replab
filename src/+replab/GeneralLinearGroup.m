classdef GeneralLinearGroup < replab.Group & replab.domain.VectorSpace
% Describes the group of square invertible real or complex matrices

    properties
        n % integer: size
        sparse % (logical): Whether to preserve sparse matrices
     end

    properties (Access = protected)
        parent % replab.domain.Matrices: General, not necessarily invertible matrices
    end

    methods

        function self = GeneralLinearGroup(field, n, sparse)
        % Constructs the group of general linear matrices
        %
        % Args:
        %   field ('R' or 'C'): Real or complex matrices
        %   n (integer): Size of the square matrices
        %   sparse (logical, optional): Whether to preserve sparse matrices, default value false
            self.field = field;
            self.n = n;
            self.parent = replab.domain.Matrices(field, n, n);
            if nargin < 3
                sparse = false;
            end
            self.sparse = sparse;
            if sparse
                self.identity = [speye(n) speye(n)];
            else
                self.identity = [eye(n) eye(n)];
            end
        end

        % Str

        function s = headerStr(self)
            s = sprintf('%d x %d invertible matrices in %s', self.n, self.n, self.field);
        end

        % Domain

        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end

        function X = sample(self)
            X = self.parent.sample;
            % a generic gaussian matrix is almost always invertible
        end

        % Semigroup/monoid/group

        function Z = compose(self, X, Y)
            Z = X * Y;
        end

        function XInv = inverse(self, X)
            XInv = inv(X);
        end

    end

end
