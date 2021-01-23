classdef UnitaryGroup < replab.CompactGroup
% Describes the group of n x n unitary (complex) matrices

    properties
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
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            realPart = randn(self.n, self.n);
            imagPart = randn(self.n, self.n);
            X = (realPart + 1i * imagPart)/sqrt(2);
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

    end

    methods % Representations

        function rep = definingRep(self)
        % Returns the defining representation of this unitary group
        %
        % The defining representation of $U(d)$ corresponds to a ``d x d`` unitary matrix.
        %
        % Returns:
        %   `.Rep`: Defining representation
            rep = replab.rep.DefiningRep(self);
        end

    end

end
