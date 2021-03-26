classdef TorusGroup < replab.CompactGroup
% Describes a compact, connected, abelian Lie group
%
% Its elements are column vectors in the space ``[0,1[^n``, i.e. floating point numbers "modulo 1", where ``n`` is the
% group dimension. The group binary operation is addition modulo 1.
%
% The elements of this torus may obey specified relations; let ``g = [g(1); g(2); ..., g(n)]`` be a group element.
% Then we would have ``mod(relations*g, 1) == 0``.
%
% There is a natural isomorphism with the group ``U(1)^n``, with image ``u`` of coefficients ``u(i) = exp(2i*pi*g(i))``.

    properties (SetAccess = protected)
        n % (integer): Torus dimension
        r % (integer): Torus rank
        relations % (integer(\*,n)): Relations obeyed by the elements of this torus
        injection % (integer(n,r)): Injection map from maximal torus to this torus
        projection % (integer(r,n)): Projection map from this torus to a maximal torus
    end

    methods

        function self = TorusGroup(relations)
            n = size(relations, 2);
            if size(relations, 1) == 0
                r = n;
                injection = eye(n);
                projection = eye(n);
            else
                [H, U] = replab.numerical.integer.hermiteNormalForm(relations');
                H = H';
                U = U';
                V = round(inv(U));
                assert(all(all(U*V == eye(n))));
                mask = all(H == 0, 1);
                injection = U(:,mask);
                projection = V(mask,:);
                r = sum(mask);
            end
            self.n = n;
            self.r = r;
            self.injection = injection;
            self.projection = projection;
            self.relations = relations;
            self.identity = zeros(n, 1);
        end

    end

    methods % Group constructions

        function S = semidirectProductFromRep(self, rep)
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
            assert(isa(H, 'replab.FiniteGroup'));
            n = rep.dimension;
            assert(n == self.n);
            for i = 1:H.nGenerators
                h = H.generator(i);
                img = rep.image(h);
                assert(all(all(round(img) == img)), 'The representation must have integer coefficients');
                assert(all(all(self.relations*img*self.injection == 0)), 'The representation must be compatible with the relations');
            end
            S = H.semidirectProduct(self, @(g, p) mod(rep.image(g)*p, 1));
        end

    end

    methods % Morphisms

        function mu = torusMorphism(self, target, torusMap)
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
            assert(target.n == n1);
            % double all(all(...)) to cater for the case of empty dimensions
            assert(all(all(target.relations * torusMap * self.injection == 0)), 'The given map is not compatible with the target relations');
            mu = self.morphismByFunction(target, @(t) mod(torusMap*t, 1), target.projection * torusMap * self.injection);
        end


        function mu = splitMorphism(self, varargin)
        % Constructs a morphism from this torus group to a direct product of torus groups
        %
        % Can be constructed either using ``blocks`` of indices, or by ``maps``; one of the arguments must be given.
        %
        % Keyword Args:
        %   blocks (cell(1,\*) of integer(1,\*)): An array of subset of indices describing which elements are mapped to which factor
        %   maps (cell(1,\*) of integer(\*,\*)): An array of torus maps describing a morphism between this torus and each of the factors
            args = struct('blocks', [], 'maps', []);
            args = replab.util.populateStruct(args, varargin);
            if iscell(args.blocks)
                assert(~iscell(args.maps), 'Only one of the ''blocks'' and ''maps'' keyword arguments must be provided');
                maps = cellfun(@(b) full(sparse(1:length(b), b, ones(1, length(b)), length(b), self.n)), args.blocks, 'uniform', 0);
            else
                assert(iscell(args.maps), 'One of ''blocks'' and ''maps'' keyword arguments must be provided');
                maps = args.maps;
            end
            dims = cellfun(@(m) size(m, 1), maps);
            torusMap = vertcat(maps{:});
            G = replab.DirectProductGroup.make(arrayfun(@(d) replab.StandardTorusGroup(d), dims, 'uniform', 0));
            mu = self.morphismByFunction(G, @(t) cellfun(@(m) mod(m*t, 1), maps, 'uniform', 0), torusMap);
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
            if ~isa(lhs, 'replab.TorusGroup') || ~isa(rhs, 'replab.TorusGroup') || lhs.n ~= rhs.n || lhs.r ~= rhs.r
                b = false;
                return
            end
            if all(all(lhs.injection == rhs.injection))
                b = true;
                return
            end
            % we need to perform column-style HNF to compare the lattices
            H1 = replab.numerical.integer.hermiteNormalForm(lhs.injection);
            H2 = replab.numerical.integer.hermiteNormalForm(rhs.injection);
            b = all(all(H1 == H2));
        end

        % Domain

        function b = eqv(self, x, y)
            b = all(x == y);
        end

        function x = sample(self)
            y = randi([0, 2^50-1], self.r, 1)/(2^50);
            x = mod(self.injection * y, 1);
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
            if self.n == self.r
                mu = replab.Isomorphism.identity(self);
            else
                T = replab.T(self.r);
                mu = T.torusMorphism(self, self.injection);
            end
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

    end

end
