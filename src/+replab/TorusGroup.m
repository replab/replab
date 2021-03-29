classdef TorusGroup < replab.CompactGroup
% Describes a compact, connected, abelian Lie group
%
% Its elements are column vectors of angles expressed in turns (see `<https://en.wikipedia.org/wiki/Turn_(angle)>_`),
% i.e. real numbers between 0 and 1, 1 not included. The group binary operation is addition modulo 1, the identity
% is the turn ``0``, and the inverse operation maps any non-zero turn ``x`` to ``1-x``.
%
% Let the coefficients ``[t1; t2; ...; tn]`` be the coefficients of this column vector, where ``n`` is the torus group
% dimension. The coefficients may or may not obey equations of the form ``mod(a1*t1 + a2*t2 + ... + an*tn, 1) == 0``,
% where the coefficients ``a1 ... an`` are integers. Those equations are collected as rows in the `.E` integer matrix.
%The size of the null space of the `.E` matrix is the number of degrees of freedom of the torus group or its rank `.r`.
%
% Such a group is isomorphic to the group of column vectors whose coefficients are unit complex numbers, where
% the identity is the number ``1``, the group binary operation is the multiplication of complex numbers and
% inverse operation the conjugation. Such vectors are written ``[u1; u2; ...; un]`` where ``ui = exp(2i*pi*ti)``.
% The equations are written as ``u1^a1 * u2^a2 * ... * un^an == 1`` with the same integer coefficients.
%
% The defining representation of this group (`.definingRep`) is the ``n x n`` matrix with diagonal elements ``u1`` ... ``un``.
%
% The explicit group elements (for example, the value returned by `.sample`) are written using the additive convention;
% while the equations are written multiplicatively in symbolic form.

    properties (SetAccess = protected)
        names % (cell(1,\*) of charstring): Names of coefficients in the multiplicative convention
        n % (integer): Torus dimension
        r % (integer): Torus rank
        E % (integer(\*,n)): Equations obeyed by the elements of this torus
        injection % (integer(n,r)): Injection map from maximal torus to this torus
        projection % (integer(r,n)): Projection map from this torus to a maximal torus
    end

    methods (Static, Access = protected)

        function names = standardNames(n)
            persistent cache
            m = max(n, 10);
            if isempty(cache)
                cache = arrayfun(@(i) sprintf('u%d', i), 1:m, 'uniform', 0);
            else
                cache = horzcat(cache, arrayfun(@(i) sprintf('u%d', i), length(cache)+1:m, 'uniform', 0));
            end
            names = cache(1:n);
        end

    end

    methods

        function self = TorusGroup(E, names)
            if nargin < 2 || isempty(names)
                n = size(E, 2);
                names = replab.TorusGroup.standardNames(n);
            else
                n = length(names);
                if isempty(E)
                    E = zeros(0, n);
                end
            end
            if size(E, 1) == 0
                r = n;
                injection = eye(n);
                projection = eye(n);
            else
                [H, U] = replab.numerical.integer.hermiteNormalForm(E');
                H = H';
                U = U';
                V = round(inv(U));
                assert(all(all(U*V == eye(n))));
                mask = all(H == 0, 1);
                injection = U(:,mask);
                projection = V(mask,:);
                r = sum(mask);
            end
            self.names = names;
            self.n = n;
            self.r = r;
            self.injection = injection;
            self.projection = projection;
            self.E = E;
            self.identity = zeros(n, 1);
        end

        function n = nEquations(self)
            n = size(self.E, 1);
        end

        function s = equation(self, i)
            e = self.E(i,:);
            l = [];
            for i = 1:length(e)
                if e(i) > 0
                    l = [l ones(1, e(i))*i];
                elseif e(i) < 0
                    l = [l -ones(1, -e(i))*i];
                end
            end
            s = replab.fp.Letters.print(l, self.names);
        end

        function T1 = subgroup(self, newE)
        % Returns a subgroup of this torus that obeys the given equations
        %
        % Args:
        %   E (integer(\*, n)): New equations
        %
        % Returns:
        %   `.TorusGroup`: Torus group satisfying the equations of this group and the ones given
            assert(size(newE, 2) == self.n);
            T1 = replab.TorusGroup([self.E; newE], self.names);
        end

        function T1 = subgroupWith(self, varargin)
            E = zeros(0, self.n);
            for i = 1:length(varargin)
                v = varargin{i};
                if isa(v, 'double')
                    E = [E; v];
                elseif ischar(v)
                    [ok, tokens] = replab.fp.Parser.lex(v, self.names);
                    assert(ok, 'Could not parse %s', v);
                    [pos, equations] = replab.fp.Parser.equation(tokens, 1);
                    if pos == 0 || tokens(1, pos) ~= replab.fp.Parser.types.END
                        error('Could not parse %s', v);
                    end
                    for i = 1:length(equations)
                        e = equations{i};
                        row = zeros(1, self.n);
                        for j = 1:length(e)
                            row(abs(e(j))) = row(abs(e(j))) + sign(e(j));
                        end
                        E = [E; row];
                    end
                else
                    error('Unsupported argument type %s', class(v));
                end
            end
            T1 = self.subgroup(E);
        end

        function T1 = withNames(self, newNames)
        % Returns a copy of this torus group with the names of the unit complex coefficients updated
        %
        % Args:
        %   newNames (cell(1,n) of charstring): Names of the complex coefficients
        %
        % Returns:
        %   `.TorusGroup`: Updated group
            assert(length(newNames) == self.n);
            T1 = replab.TorusGroup(self.E, newNames);
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
                assert(all(all(self.E*img*self.injection == 0)), 'The representation must be compatible with the group equations');
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
            assert(all(all(target.E * torusMap * self.injection == 0)), 'The given map is not compatible with the target equations');
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
            G = replab.DirectProductGroup.make(arrayfun(@(d) replab.T(d), dims, 'uniform', 0));
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
            rep = self.diagonalRep(eye(self.n));
        end

        function rep = diagonalRepWith(self, varargin)
        % Returns a diagonal rep of this torus group
            torusMap = zeros(0, self.n);
            for i = 1:length(varargin)
                v = varargin{i};
                if ~ischar(v)
                    error('Unsupported argument type %s', class(v));
                end
                assert(ischar(v));
                [ok, tokens] = replab.fp.Parser.lex(v, self.names);
                assert(ok, 'Could not parse %s', v);
                [pos, res] = replab.fp.Parser.word(tokens, 1);
                if pos == 0 || tokens(1, pos) ~= replab.fp.Parser.types.END
                    error('Could not parse %s', v);
                end
                row = zeros(1, self.n);
                for j = 1:length(res)
                    row(abs(res(j))) = row(abs(res(j))) + sign(res(j));
                end
                torusMap = [torusMap; row];
            end
            rep = self.diagonalRep(torusMap);
        end

        function rep = diagonalRep(self, torusMap)
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

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.CompactGroup(self);
            names{1,end+1} = 'E';
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.CompactGroup(self);
            for i = 1:self.nEquations
                names{1, end+1} = sprintf('equation(%d)', i);
                values{1, end+1} = self.equation(i);
            end
        end

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
