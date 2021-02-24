classdef TorusGroup < replab.CompactGroup
% Describes a compact, connected, abelian Lie group
%
% Its elements are row vectors in the space ``[0,1[^n``, i.e. floating point numbers "modulo 1"; the group binary operation
% is addition.
%
% Though it is not enforced, those numbers should be of the form ``m * 2^-50`` where ``m`` is a 50-bit unsigned integer, so
% that composition can be computed
%
% Let ``g = [g(1), g(2), ..., g(n)]`` be a group element. There is an isomorphism with the group ``U(1)^n``, whose image ``u``
% has coefficients: ``u(i) = exp(2i*pi*g(i))``.

    properties (SetAccess = protected)
        n % (integer): Torus dimension
    end

    methods

        function self = TorusGroup(n)
            self.n = n;
            self.identity = zeros(1, n);
        end

        function X = toMatrix(self, x)
            X = diag(exp(2i*pi*x));
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(x == y);
        end

        function x = sample(self)
            x = randi([0, 2^50-1], 1, self.n)/(2^50);
        end

        % Monoid

        function z = compose(self, x, y)
            z = mod(x + y, 1);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = mod(-x, 1);
        end

    end

end
