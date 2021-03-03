classdef TorusGroup < replab.CompactGroup
% Describes a compact, connected, abelian Lie group
%
% Its elements are column vectors in the space ``[0,1[^n``, i.e. floating point numbers "modulo 1"; the group binary operation
% is addition.
%
% Though it is not enforced, those numbers should be of the form ``m * 2^-50`` where ``m`` is a 50-bit unsigned integer,
% so that composition can be computed exactly.
%
% Let ``g = [g(1); g(2); ..., g(n)]`` be a group element.
% There is an isomorphism with the group ``U(1)^n``, whose image ``u`` has coefficients: ``u(i) = exp(2i*pi*g(i))``.

    properties (SetAccess = protected)
        n % (integer): Torus dimension
    end

    methods

        function self = TorusGroup(n)
            self.n = n;
            self.identity = zeros(n, 1);
        end

        function X = toMatrix(self, x)
            assert(length(x) == self.n);
            X = diag(exp(2i*pi*x));
        end

        function rep = definingRep(self)
        % Returns the defining representation of this torus group
        %
        % The defining representation of a torus group of dimension ``n`` is a n-dimensional representation,
        % whose images are diagonal in the standard Euclidean basis.
        %
        % Returns:
        %   `.Rep`: Representation
            rep = self.mapRep(eye(self.n));
        end

        function rep = mapRep(self, torusMap)
        % Returns a diagonal representation of this torus group
        %
        % Constructs a representation using a ``torusMap``, whose images are diagonal matrices. The diagonal
        % elements are monomials of the phases of the torus group. Each row of ``torusMap`` corresponds to
        % a diagonal element, and the columns give the integer exponents corresponding to every phase.
        %
        % Args:
        %   torusMap (integer(d, n)): Phase exponents
        %
        % Returns:
        %   `.Rep`: Representation
            rep = replab.rep.TorusRep(self, torusMap);
        end

    end

    methods % Implementations

        % Obj

        function b = eq(lhs, rhs)
            b = isa(lhs, 'replab.TorusGroup') && isa(rhs, 'replab.TorusGroup') && lhs.n == rhs.n;
        end

        % Domain

        function b = eqv(self, x, y)
            b = all(x == y);
        end

        function x = sample(self)
            x = randi([0, 2^50-1], self.n, 1)/(2^50);
        end

        % Monoid

        function z = compose(self, x, y)
            z = mod(x + y, 1);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = mod(-x, 1);
        end

        % CompactGroup

        function b = hasReconstruction(self)
            b = true;
        end

        function [mu, R] = reconstruction(self)
            R = replab.SetProduct.identity(self);
            mu = replab.Isomorphism.identity(self);
        end

    end

end
