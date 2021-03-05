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

    end

    methods % Morphisms

        function mu = torusMorphism(self, torusMap)
        % Constructs a morphism from this torus group to another torus group from an integer matrix
        %
        % Args:
        %   torusMap (integer(n1,n)): Integer matrix describing the action on the torus elements
        %
        % Returns:
        %   `.Morphism`: A morphism from this torus group to a torus group of dimension ``n1``
            n = size(torusMap, 2);
            n1 = size(torusMap, 1);
            assert(self.n == n);
            T1 = replab.TorusGroup(n1);
            mu = self.morphismByFunction(T1, @(t) torusMap*t, torusMap);
        end

    end

    methods % Representations

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

    methods (Static)

        function rho = torusRepRealImage(t)
        % Returns the sparse block-diagonal matrix corresponding to the given torus element
        %
        % The element ``t(i)`` will correspond to a 2x2 block: ``[c s; -s c]`` where
        % ``c = cos(2*pi*t(i))`` and ``s = sin(2*pi*t(i))``.
        %
        % Args:
        %   t (double(n,1)): Torus group element
        %
        % Returns:
        %   double(2*n,2*n): Block diagonal matrix of 2x2 real matrices
            n = length(t);
            blocks = arrayfun(@(ti) sparse([cos(2*pi*ti) -sin(2*pi*ti); sin(2*pi*ti) cos(2*pi*ti)]), t, 'uniform', false);
            rho = blkdiag(blocks{:});
        end

        function rho = torusRepImage(t)
        % Returns the sparse diagonal matrix corresponding to the given torus element
        %
        % Equivalent to ``torusGroup.definingRep.image(t)``, without constructing additional objects.
        %
        % Args:
        %   t (double(n,1)): Torus group element
        %
        % Returns:
        %   double(n,n), sparse: Diagonal matrix corresponding to the torus element
            n = length(t);
            rho = sparse(1:n, 1:n, exp(2i*pi*mod(t', 1)));
        end

        function S = semidirectProductFromRep(rep)
        % Returns the semidirect of a finite group on the torus elements using one of its integer representations
        %
        % Invertible integer matrices represent automorphisms of a torus group of the same dimension; thus, given a
        % group ``G``, a representation of ``G`` with integer images defines a morphism from ``G`` to the automorphism
        % group of a torus ``T``.
        %
        % One can then use this morphism to construct an outer semidirect product of ``G`` acting on ``T``.
        %
        % Args:
        %   rep (`.Rep`): Group representation with integer coefficients
        %
        % Returns:
        %   `.SemidirectProduct`: Semidirect product of ``rep.group`` and the torus group
            H = rep.group;
            n = rep.dimension;
            assert(isa(H, 'replab.FiniteGroup'));
            for i = 1:H.nGenerators
                h = H.generator(i);
                img = rep.image(h);
                assert(all(all(round(img) == img)), 'The representation must have integer coefficients');
            end
            T = replab.TorusGroup(n);
            S = H.semidirectProduct(T, @(g, p) mod(rep.image(g)*p, 1));
        end

    end

end
