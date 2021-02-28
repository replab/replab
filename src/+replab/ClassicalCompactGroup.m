classdef ClassicalCompactGroup < replab.CompactGroup
% Describes a classical compact group of square real/complex/quaternionic matrices
%
% Group elements are matrices whose coefficients are encoded using double floating-point arithmetic. Of course, they will only be
% approximately orthogonal and thus do not represent directly a group element.
%
% Instead, we postulate that the matrix encodes the nearest group element in Frobenius norm. For non-special groups, this element
% can be found as ``U = X*sqrtm(X'*X)`` where ``X`` is the approximate element.
%
% While inverses should be exact (we take the conjugate transpose), the composition is only approximate. Thus, to be able to check
% the usual group laws, we compare group elements up to a certain tolerance, defined in `+replab.+globals.matrixGroupTol`.

    properties (SetAccess = protected)
        n % (integer): Dimension of the matrices
        algebra % ('R', 'C' or 'H'): Unital associative algebra over which the group is defined
        isSpecial % (logical): Whether the matrices should have determinant one; must be false when algebra = 'H'
    end

    methods

        function self = ClassicalCompactGroup(n, algebra, isSpecial)
        % Constructs a classical compact group
        %
        % Args:
        %   n (integer): Size of the ``n x n`` matrices
        %   algebra ('R', 'C' or 'H'): Unital associatvie algebra over which the group is defined
        %   isSpecial (logical): Whether the matrices should have ``det(X) == 1``
            assert(ischar(algebra) && length(algebra) == 1 && any(algebra == 'RCH'));
            assert(algebra ~= 'H' || ~isSpecial);
            self.n = n;
            self.algebra = algebra;
            self.isSpecial = logical(isSpecial);
            if algebra == 'H'
                self.identity = replab.H(eye(n));
            else
                self.identity = eye(n);
            end
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            switch self.algebra
              case 'R'
                if self.isSpecial
                    s = sprintf('Special orthogonal group SO(%d)', self.n);
                else
                    s = sprintf('Orthogonal group O(%d)', self.n);
                end
              case 'C'
                if self.isSpecial
                    s = sprintf('Special unitary group SU(%d)', self.n);
                else
                    s = sprintf('Unitary group U(%d)', self.n);
                end
              case 'H'
                s = sprintf('Compact symplectic group Sp(%d)', self.n);
            end
        end

        % Domain

        function b = eqv(self, X, Y)
            b = norm(X - Y, 'fro') <= replab.globals.matrixGroupTol;
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
            [X, detX] = replab.numerical.randomUnitaryOver(self.n, self.algebra);
            if self.isSpecial
                detX = sign(detX);
                c = randi(self.n);
                X(:,c) = X(:,c)/detX;
            end
        end

        function b = hasReconstruction(self)
            if any(self.algebra == 'CH') && ~self.isSpecial
                b = true;
            else
                b = false;
            end
        end

        function [mu, R] = reconstruction(self)
            assert(self.hasReconstruction);
            R = replab.SetProduct.identity(self);
            T = replab.TorusGroup(self.n);
            if self.algebra == 'C'
                mu = T.morphismByFunction(self, @(x) T.toMatrix(x));
            else
                mu = T.morphismByFunction(self, @(x) replab.H(T.toMatrix(x)));
            end
        end

    end

    methods % Representations

        function rep = definingRep(self, field)
        % Returns the defining representation of this unitary group
        %
        % The defining representation of classical compact group over 'R' or 'C' is a real/complex
        % real/complex ``n x n`` orthogonal or unitary matrix. The defining representation of a
        % classical group over 'H' corresponds to a complex matrix encoding the quaternion division
        % algebra.
        %
        % Args:
        %   field ('R', 'C', optional): Field over which to define the matrices, default: C
        %
        % Returns:
        %   `.Rep`: Defining representation
            if nargin < 2 || isempty(field)
                if self.algebra == 'H'
                    field = 'C';
                else
                    field = self.algebra;
                end
            end
            rep = replab.rep.DefiningRep(self, field);
        end

    end

end
