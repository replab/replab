classdef LeftCosets < replab.Cosets
% Describes the set of left cosets of a finite group
%
% Let $H$ be a subgroup of a group $G$. Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% The set of such left cosets is often written $G / H = \{ g H : g \in G \}$.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.LeftCosetsLaws(self);
        end

    end

    methods

        function t = cosetRepresentative(self, g)
        % Returns the canonical coset representative corresponding to the given element
        %
        % If ``t = cosetRepresentative(g)``, then $t^{-1} g = h \in H$, with decomposition $g = t h$.
        %
        % Moreover ``L.cosetRepresentative(g) == L.cosetRepresentative(compose(g, h))`` for any $h \in H$.
        %
        % Finally, ``L.cosetRepresentative(g) == L.coset(g).representative``.
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   t (element of `.group`): Coset canonical representative
            error('Abstract');
        end

        function C = elements(self)
        % Returns the set of left cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.LeftCoset`: Set of left cosets
            C = cellfun(@(t) self.subgroup.leftCoset(t, 'group', self.group), self.transversal, 'uniform', 0);
        end

        function mu = leftAction(self)
        % Returns, as a morphism, the action of the given group of its left cosets
        %
        % Returns:
        %   `.FiniteMorphism`: Morphism of this group into a permutation group
            error('Abstract');
        end

        function s = nElements(self)
        % Returns the number of left cosets
        %
        % Returns:
        %   integer: Number of left cosets
            s = self.group.order / self.subgroup.order;
            assert(s <= 2^53 - 1);
            s = double(s);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            error('Abstract');
        end

    end

end
