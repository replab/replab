classdef TorusGroup < replab.CompactGroup
% Describes a compact, connected, abelian Lie group
%
% Its elements are column vectors of angles expressed in turns (see `<https://en.wikipedia.org/wiki/Turn_(angle)>_`),
% i.e. real numbers between 0 and 1, 1 not included. The group binary operation is addition modulo 1, the identity
% is the turn ``0``, and the inverse operation maps any non-zero turn ``x`` to ``1-x``.
%
% When there are no additional restrictions on the elements of a given torus group, it is the standard torus group,
% and it can be constructed with the function `+replab.T`.
%
% Torus morphisms from integer matrices: `<https://math.stackexchange.com/questions/2229088/automorphisms-of-a-torus>_`.
% Let ``s = [s1;...;sn]`` be the element of a (standard) torus of dimension ``n`` and ``t = [t1;...;tm]`` the element of
% a (standard) torus of dimension ``m`. A torus morphism ``f`` is given by a ``m x n`` integer matrix ``M``, with the
% image of ``s`` being written ``t = f(s) = mod(M * s, 1)``.
%
% Let the coefficients ``[t1; t2; ...; tn]`` be the coefficients of this column vector, where ``n`` is the torus group
% dimension. The coefficients may or may not obey equations of the form ``mod(a1*t1 + a2*t2 + ... + an*tn, 1) == 0``,
% where the coefficients ``a1 ... an`` are integers. Those equations are collected as rows in the `.E` integer matrix.
% Equivalently, the torus group obeying those relations is the kernel of the torus morphism encoded in `.E`.
%
% How do we deal with those equations? We use the Hermite normal form of the matrix ``E`` to extract a pair of morphisms
% from ``T`` (whose elements obey the equations) to a standard torus ``S`` (whose dimension is then the rank of ``T``).
% The morphism `.projection` from ``T`` to ``S`` is surjective, and the morphism `.injection` from ``S`` to ``T`` is injective,
% and while ``projection(injection(s)) = s``. Both morphisms are given by integer matrices.
%
% The defining representation of this group (`.definingRep`) is the ``n x n`` matrix with diagonal elements ``u1`` ... ``un``.
%
% Automorphisms of a torus group are given by invertible square integer matrices ``M`` such that ``E*M*injection = 0``.
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
        % Constructs a torus group obeying the given equations, with optional element names
        %
        % Args:
        %   E (integer(\*,n)): Integer matrix describing the equations; its number of columns provide the torus dimension
        %   names (cell(1,n) of charstring, optional): Names of the group element coefficients in their unit complex form
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
        % Returns the number of equations defining this torus group
        %
        % Returns:
        %   integer: Number of equations
            n = size(self.E, 1);
        end

        function s = equation(self, i)
        % Returns the i-th equation defining this torus group in monomial form
        %
        % Args:
        %   i (integer): Equation index
        %
        % Returns:
        %   charstring: Equation in monomial form (whose value must be equal to one)
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

        function s = equations(self)
        % Returns all equations that this torus group obeys in monomial form
        %
        % Returns:
        %   cell(1,\*) of charstring: Equations defining this torus group in monomial form
            s = arrayfun(@(i) self.equation(i), 1:self.nEquations, 'uniform', 0);
        end

        function T1 = subgroup(self, newE)
        % Returns a subgroup of this torus that obeys the equations given in integer matrix form
        %
        % Args:
        %   E (integer(\*, n)): New equations
        %
        % Returns:
        %   `.TorusGroup`: Torus group satisfying the equations of this group and the additional ones given
            assert(size(newE, 2) == self.n);
            T1 = replab.TorusGroup([self.E; newE], self.names);
        end

        function T1 = subgroupWith(self, varargin)
        % Returns a subgroup of this torus that obeys the equations given in monomial form
        %
        % This function accepts a variable number of arguments. Each argument must be of the form:
        %
        % * An integer matrix with `.n` columns, defining equations as in the call to `.subgroup`,
        % * A charstring containing an equation in monomial form. If the equation has no right-hand side, the
        %   right-hand side is assumed to be the identity.
        %
        % Returns:
        %   `.TorusGroup`: Torus group satisfying the equations of this group and the additional ones given
            E = self.parseTorusMap(varargin, 'equations');
            T1 = self.subgroup(E);
        end

        function T1 = withNames(self, newNames)
        % Returns a copy of this torus group with the names of the unit complex coefficients replaced
        %
        % Args:
        %   newNames (cell(1,n) of charstring): New names of the complex coefficients
        %
        % Returns:
        %   `.TorusGroup`: Updated group
            assert(length(newNames) == self.n);
            T1 = replab.TorusGroup(self.E, newNames);
        end

    end

    methods % Group constructions

        function S = semidirectProductByFiniteGroup(self, group, varargin)
        % Returns the semidirect product of a finite group and this torus, with the finite group action defined by images
        %
        % Args:
        %   group (`.FiniteGroup`): Finite group acting on this torus group
        %
        % Keyword Args:
        %   preimages (cell(1, \*) of ``group`` elements): Preimages of the morphism which generate ``group``, defaults to ``group.generators``
        %   images (cell(1, \*) of ``target`` elements): Action on the torus for the given preimages, defaults to a permutation action if ``group`` is a permutation group
        %
        % Returns:
        %   `.SemidirectProduct`: Semidirect product of ``group`` and this torus group
            assert(isa(group, 'replab.FiniteGroup'));
            args = struct('preimages', {group.generators}, 'images', {[]});
            args = replab.util.populateStruct(args, varargin);
            n = self.n;
            r = self.r;
            if group.isTrivial
                rep = group.trivialRep('R', n);
            elseif isempty(args.images)
                assert(isa(group, 'replab.PermutationGroup'), 'If the images argument is missing, the group must be a permutation group');
                rep = group.naturalRep;
            else
                if isempty(args.preimages)
                    args.preimages = group.generators;
                end
                N = length(args.preimages);
                assert(length(args.images) == N, 'As many images as preimages must be given');
                repImages = cell(1, N);
                for i = 1:N
                    m = args.images{i};
                    assert(isa(m, 'replab.Isomorphism') && m.source == self && m.target == self, 'Images must be automorphisms of this torus group');
                    repImages{i} = self.injection*m.torusMap*self.projection;
                end
                rep = group.repByImages('R', n, 'preimages', args.preimages, 'images', repImages);
            end
            S = self.semidirectProductFromRep(rep);
        end

        function S = semidirectProductFromRep(self, rep)
        % Returns the semidirect product of a finite group on the torus elements using one of its integer representations
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
        %   `.SemidirectProduct`: Semidirect product of ``rep.group`` and this torus group
            H = rep.group;
            assert(isa(H, 'replab.FiniteGroup'));
            assert(isa(rep, 'replab.Rep') && rep.overR);
            n = rep.dimension;
            assert(n == self.n, 'Dimension of the action not matching the dimension of the torus');
            for i = 1:H.nGenerators
                h = H.generator(i);
                img = rep.image(h);
                invImg = rep.inverseImage(h);
                assert(all(all(round(img) == img)), 'The representation must have integer coefficients');
                assert(all(all(round(invImg) == invImg)), 'The representation must have integer coefficients');
                assert(all(all(self.E*img*self.injection == 0)), 'The representation must be compatible with the group equations');
                assert(all(all(self.E*invImg*self.injection == 0)), 'The representation must be compatible with the group equations');
            end
            S = H.semidirectProduct(self, @(g, p) mod(rep.image(g)*p, 1));
        end

    end

    methods (Access = protected)

        function M = parseTorusMap(self, args, type)
        % Parses a series of equations or monomials
        %
        % Args:
        %   args (cell(1,\*) of either integer(\*,n) or charstring): Arguments to parse
        %   type ('equations', 'monomials'): Whether to allow equations (such as 'x*y=y*z=1')
            M = zeros(0, self.n);
            for i = 1:length(args)
                a = args{i};
                if isa(a, 'double')
                    assert(size(a, 2) == self.n);
                    assert(all(all(round(a) == a)), 'Matrix must be integer-valued');
                    M = [M; a];
                elseif ischar(a)
                    [ok, tokens] = replab.fp.Parser.lex(a, self.names);
                    assert(ok, 'Could not parse %s', a);
                    switch type
                      case 'equations'
                        [pos, words] = replab.fp.Parser.equation(tokens, 1);
                      case 'monomials'
                        [pos, word] = replab.fp.Parser.word(tokens, 1);
                        words = {word};
                      otherwise
                        error('Type %s not supported', type);
                    end
                    if pos == 0 || tokens(1, pos) ~= replab.fp.Parser.types.END
                        error('Could not parse %s', a);
                    end
                    for i = 1:length(words)
                        w = words{i};
                        row = zeros(1, self.n);
                        for j = 1:length(w)
                            row(abs(w(j))) = row(abs(w(j))) + sign(w(j));
                        end
                        M = [M; row];
                    end
                else
                    error('Unsupported argument type %s', class(v));
                end
            end
        end

    end

    methods % Morphisms

        function mu = automorphism(self, varargin)
        % Constructs an automorphism of this torus from an integer matrix describing a torus map
        %
        % This function accepts a variable number of arguments. Each argument must be of the form:
        %
        % * An integer matrix with `.n` columns, defining equations as in the call to `.subgroup`,
        % * A charstring containing an equation in monomial form.
        %
        % Those arguments will be parsed and concatenated vertically to form an integer matrix
        % describing the torus map. That matrix must be invertible.
        %
        % Returns:
        %   `.Isomorphism`: The constructed automorphism
            M = self.parseTorusMap(varargin, 'monomials');
            assert(size(M, 1) == self.n);
            assert(abs(det(M)) == 1, 'Matrix must be invertible for an automorphism');
            invM = round(inv(M));
            assert(all(all(self.E * M * self.injection == 0)), 'Map must respect the torus structure');
            assert(all(all(self.E * invM * self.injection == 0)), 'Map must respect the torus structure');
            mu = replab.Isomorphism.lambda(self, self, @(x) mod(invM*x, 1), @(x) mod(M*x, 1), self.projection*M*self.injection);
        end

        function mu = endomorphism(self, varargin)
        % Constructs an endomorphism of this torus from an integer matrix describing a torus map
        %
        % This function accepts a variable number of arguments. Each argument must be of the form:
        %
        % * An integer matrix with `.n` columns, defining equations as in the call to `.subgroup`,
        % * A charstring containing an equation in monomial form.
        %
        % Those arguments will be parsed and concatenated vertically to form an integer matrix
        % describing the torus map.
        %
        % Returns:
        %   `.Morphism`: The constructed endomorphism
            M = self.parseTorusMap(varargin, 'monomials');
            assert(all(all(self.E * M * self.injection == 0)), 'Map must respect the torus structure');
            assert(size(M, 1) == self.n);
            mu = replab.Morphism.lambda(self, self, @(x) mod(M*x, 1), self.projection*M*self.injection);
        end

        function mu = torusMorphism(self, target, varargin)
        % Constructs a morphism from this torus group to another torus group from an integer matrix
        %
        % After the ``target`` argument, this function accepts a variable number of arguments.
        % Each of these arguments must be of the form:
        %
        % * An integer matrix with `.n` columns, defining equations as in the call to `.subgroup`,
        % * A charstring containing an equation in monomial form.
        %
        % Args:
        %   target (`.TorusGroup`): Torus group target of the morphism
        %
        % Returns:
        %   `.Morphism`: A morphism from this torus group to a torus group of dimension ``n1``
            torusMap = self.parseTorusMap(varargin, 'monomials');
            n = size(torusMap, 2);
            n1 = size(torusMap, 1);
            assert(self.n == n);
            assert(target.n == n1);
            % double all(all(...)) to cater for the case of empty dimensions
            assert(all(all(target.E * torusMap * self.injection == 0)), 'The given map is not compatible with the target equations');
            mu = self.morphismByFunction(target, @(t) mod(torusMap*t, 1), target.projection * torusMap * self.injection);
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

        function h = headerStr(self)
            h = sprintf('Torus group of dimension n=%d and rank r=%d', self.n, self.r);
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.CompactGroup(self);
            names{1,end+1} = 'E';
            names{1,end+1} = 'n';
            names{1,end+1} = 'r';
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
