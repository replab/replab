classdef SemidirectProductGroup < replab.Group
% Semidirect product of compact groups through the external/outer construction
%
% This is an abstract base class. Call `.CompactGroup.semidirectProduct` or `.make` to construct an instance.
%
% Construction and multiplication rule
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% An outer semidirect product constructs a group from two given groups, with its set of elements being the Cartesian
% product of those two groups, and the multiplication operation given by the action of one group on the other.
%
% Formally, we have two groups, $H$ and $N$, with $H$ acting on $N$ in the following way. We have a group
% action $\phi: H \times N \rightarrow N$, obeying the following:
%
% - $\phi(h_1, \phi(h_2, n)) = \phi(h_1 h_2, n)$
% - $\phi(h, n_1) \phi(h, n_2) = \phi(h, n_1 n_2)$
%
% for $h,h_1,h_2 \in H$ and $n,n_1,n_2 \in N$.
%
% Equivalently, this action corresponds to a group morphism $\phi: H \rightarrow (N \rightarrow N)$ or
% $\phi: H \rightarrow \operatorname{Aut}(N)$ from $H$ into the automorphism group of $N$.
%
% We now construct the semidirect product $G = H \ltimes N$. An element of the product is written $G = (h, n)$,
% where $h \in H$ and $n \in N$.
%
% Note that compared to other references such as Wikipedia ( `<https://en.wikipedia.org/wiki/Semidirect_product>`_ ),
% we write the element of the acting group $H$ first, and this changes the rule for multiplying group elements.
%
% The multiplication rule of the group is $(h_1, n_1) \cdot (h_2, n_2) = (h_1 h_2, \phi(h_2^{-1}, n_1) n_2)$.
% The inverse of $(h, n)$ is given by $(h^{-1}, \phi(h, n))$.
%
% How do we interpret this multiplication rule? First, we want to interpret the elements of $H$ and $N$ as "being part" of $G$.
%
% We can embed elements of $H$ into $G$ by writing the injection $f(h) = (h, 1_N)$ where $1_N$ is the identity of $N$.
% We can embed elements of $N$ into $G$ by writing the injection $g(n) = (1_H, n)$ where $1_H$ is the identity of $H$.
%
% Both $f$ and $g$ are morphisms: $f(h_1) f(h_2) = f(h_1 h_2)$ and $g(n_1) g(n_2) = g(n_1 n_2)$.
%
% Those those injections are available using the methods `.Hinjection` and `.Ninjection`.
%
% Noting that $(h,n) = f(h) g(n)$, we write:
%
% $(h_1,n_1) (h_2,n_2) = f(h_1) g(n_1) f(h_2) g(n_2) = h_1 n_1 h_2 n_2$, with the $f()$ and $g()$ injections implicit.
%
% If we could swap $n_1$ and $h_2$, as is the case with a direct product, we could rewrite the product as a pair of elements of
% $H$ and $N$ respectively. Now, we write $h_1 n_1 h_2 n_2 = h_1 h_2 h_2^{-1} n_1 h_2 n_2$, and we define
% $h_2^{-1} n_1 h_2 = \phi(h_2^{-1}, n_1)$.
%
% Thus, $h_1 n_1 h_2 n_2 = (h_1 h_2) (\phi(h_2^{-1}, n_1) n_2) = (h_1 h_2, \phi(h_2^{-1}, n_1) n_2)$.
%
% Thus, we interpret that group action $\phi$ as defining the conjugation of elements of $N$ by elements of $H$.
%
% Representations
% ~~~~~~~~~~~~~~~
%
% Representations of semidirect products can be constructed using the `.semidirectProductRep` method, given compatible
% representations of the groups $H$ and $N$.
%
% Implementation note
% ~~~~~~~~~~~~~~~~~~~
%
% As semidirect product groups are used as a base for wreath product groups,  the constructors are duplicated in subclasses
% as to keep a simple hierarchy of constructor calls.
%
% Example:
%   >>> N = replab.S(3);
%   >>> H = replab.PermutationGroup.cyclic(3);
%   >>> A = N.innerAutomorphism([2 3 1]);
%   >>> phi = H.morphismByImages(replab.AutomorphismGroup(N), 'images', {A});
%   >>> sd = H.semidirectProduct(N, @(h, n) phi.imageElement(h).imageElement(n));
%   >>> sd.laws.checkSilent
%       1

    properties (SetAccess = protected)
        H % (`+replab.CompactGroup`): Group acting
        N % (`+replab.CompactGroup`): Group acted upon
        phi % (`+replab.Action`): Action of H on N
    end

    methods (Static) % SemidirectProductGroup creation

        function prd = make(phi)
        % Constructs a semidirect product group from an action
        %
        % Args:
        %   phi (`+replab.Action`): Action of a compact group on another compact group
        %
        % Returns:
        %   `.SemidirectProductGroup`: A specialized instance of `.SemidirectProductGroup`
            isFinite = isa(phi.G, 'replab.FiniteGroup') && isa(phi.P, 'replab.FiniteGroup');
            if isFinite
                prd = replab.prods.SemidirectProductGroup_finite(phi, 'self');
            else
                prd = replab.prods.SemidirectProductGroup_compact(phi);
            end
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end

        function g = sample(self)
            g = {self.H.sample self.N.sample};
        end

        % Monoid

        function z = compose(self, x, y)
            xh = x{1};
            xn = x{2};
            yh = y{1};
            yn = y{2};
            % Relation to phi is the conjugation
            % phi_h(n) = h n h^-1
            % we have z = xh xn yh yn = xh yh yh^-1 xn yh yn =
            % = xh yh phi_(yh^-1)(xn) yn
            % and thus
            % zh = xh yh
            % zn = phi_(yh^-1)(xn) yn
            yhinv = self.H.inverse(yh);
            zh = self.H.compose(xh, yh);
            zn = self.N.compose(self.phi.leftAction(yhinv, xn), yn);
            z = {zh zn};
        end

        % Group

        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = self.N.inverse(self.phi.leftAction(xh, xn));
            z = {zh zn};
        end

    end

    methods % Morphisms

        function m = Hinjection(self)
        % Returns the morphism that expresses elements of `.H` in the semidirect product
        %
        % Returns:
        %   `.Morphism`: The morphism from `.H` to this semidirect product group
            m = self.H.morphismByFunction(self, @(g) {g, self.N.identity});
        end

        function m = Hprojection(self)
        % Returns the morphism that projects elements of the semidirect product into `.H`
        %
        % Note that for a semidirect product ``G``, the morphism ``G.Hinjection.andThen(G.Hprojection``
        % is the identity.
        %
        % Returns:
        %   `.Morphism`: The morphism from the semidirect product to `.H`
            m = self.morphismByFunction(self.H, @(g) g{1});
        end

        function m = Ninjection(self)
        % Returns the morphism that expresses elements of `.N` in the semidirect product
        %
        % Returns:
        %   `.Morphism`: The morphism from `.N` to this semidirect product group
            if self.hasReconstruction
                r = self.N.maximalTorusDimension;
                tm = eye(r);
            else
                tm = [];
            end
            m = self.N.morphismByFunction(self, @(g) {self.H.identity, g}, tm);
        end

    end

    methods % Representations

        function rep = semidirectProductRep(self, Hrep, Nrep)
        % Constructs a representation by a product of representations of the two groups in the semidirect product construction
        %
        % Those representations ``Hrep`` and ``Nrep`` must obey the following law.
        %
        % Let ``phi`` be the semidirect product group uncurried homomorphism with ``h \in H`` and ``n \in N``:
        %
        % ``Nrep.image(phi(h, n)) == Hrep.image(h) * Nrep.image(n) * Hrep.inverseImage(h)``.
        %
        % Args:
        %   Hrep (`.Rep`): Representation of `.H`
        %   Nrep (`.Rep`): Representation of `.N`
        %
        % Returns:
        %   `.Rep`: Group representation
            assert(Hrep.dimension == Nrep.dimension);
            assert(Hrep.field == Nrep.field);
            rep = replab.prods.SemidirectProductRep(self, Hrep, Nrep);
        end

    end

end
