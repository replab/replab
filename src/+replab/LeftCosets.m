classdef LeftCosets < replab.Str
% Describes the set of left right cosets of a nice finite group
%
% Let $H$ be a subgroup of a group $G$. Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% The set of such left cosets is often written $G / H = \{ g H : g \in G \}$.
%
% RepLAB works with left cosets by translating the operations to the ones on the right cosets.

    properties (SetAccess = protected)
        group % (`.NiceFiniteGroup`): Group
        subgroup % (`.NiceFiniteGroup`): Subgroup of `.group`
    end

    methods

        function t = canonicalRepresentative(self, g)
        % Returns the coset representative corresponding to the given element
        %
        % It is thus guaranteed that if ``t = canonicalRepresentative(g)``, then $t^{-1} g = h \in H$, with
        % the decomposition $g = t h$.
        %
        % Moreover ``R.canonicalRepresentative(g) == R.canonicalRepresentative(compose(g, h))`` for any $h \in H$.
        %
        % For permutation groups, the canonical representative returned is the one which is lexicographically
        % minimal in its coset.
        %
        % Args:
        %   g (element of `group`): Group element
        %
        % Returns:
        %   t (element of `group`): Coset canonical representative
            error('Abstract');
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
            error('Abstract');
        end

        function C = coset(self, g)
        % Returns all the elements in the coset containing the given element
            C = cellfun(@(h) self.group.compose(g, h), self.subgroup.elements.toCell, 'uniform', 0);
        end

        function C = cosets(self)
        % Returns the set of right cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of cell(1,\*) of elements of `.group`: Set of right cosets
            T = self.transversal;
            C = cellfun(@(t) self.coset(t), T, 'uniform', 0);
        end

        function mu = action(self)
        % Returns a morphism that describes the action of the group on cosets
            error('Abstract');
        end

    end

end
