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
            b = true;
        end

        function [mu, R] = reconstruction(self)
            n = self.n;
            if self.algebra == 'R' && ~self.isSpecial
                % O(n) has two connected components
                n1 = round((n+1)/2);
                if mod(n1, 2) == 0
                    n1 = n1 - 1;
                end
                flip = diag([ones(1, n1) -ones(1, n-n1)]);
                R = replab.SetProduct(self, {{self.identity flip}}, true);
            else
                % U(n), SU(n), SO(n), Sp(n) have one connected component
                R = replab.SetProduct.identity(self);
            end
            switch self.algebra
              case 'R'
                if mod(n, 2) == 0
                    T = replab.StandardTorusGroup(n/2);
                    mu = T.morphismByFunction(self, @(t) replab.TorusGroup.torusRepRealImage(t), eye(n/2));
                else
                    T = replab.StandardTorusGroup((n-1)/2);
                    mu = T.morphismByFunction(self, @(t) blkdiag(replab.TorusGroup.torusRepRealImage(t), sparse(1)) , eye((n-1)/2));
                end
              case 'C'
                if self.isSpecial
                    T = replab.StandardTorusGroup(n-1);
                    torusMap = [eye(n-1); -ones(1, n-1)];
                    mu = T.morphismByFunction(self, @(t) replab.TorusGroup.torusRepImage(torusMap * t), eye(n-1));
                else
                    T = replab.StandardTorusGroup(n);
                    mu = T.morphismByFunction(self, @(t) replab.TorusGroup.torusRepImage(t), eye(n));
                end
              case 'H'
                T = replab.StandardTorusGroup(n);
                mu = T.morphismByFunction(self, @(t) replab.H(replab.TorusGroup.torusRepImage(t)), eye(n));
            end
        end

    end

    methods % Representations

        function rep = definingRep(self, field)
        % Returns the defining representation of this unitary group
        %
        % "Defining representation" is a loosely defined term designing the smallest faithful representation of a
        % group. In our case, the defining representation of ``O(n)`` or ``SO(n)`` is a n-dimensional real representation,
        % the defining representation of ``U(n)`` or ``SU(n)`` is a n-dimension complex representation, while
        % the defining representation of ``Sp(n)`` is a 2n-dimension complex representation.
        %
        % All these representations are unitary. By specifying the field, one can encode the representation in a
        % different field: if the natural field is real (``O(n)``, ``SO(n)``), the representation can be complexified
        % into a representation of the same dimension. If the natural field is complex (``U(n)``, ``SU(n)``, ``Sp(n)``),
        % the complex representation can be encoded into a real representation, doubling its dimension.
        %
        % Args:
        %   field ('R', 'C', optional): Field over which to define the matrices, default: group dependent, see above
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
