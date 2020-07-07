classdef RightCosets < replab.Str
% Describes the set of right cosets of a nice finite group
%
% Let $H$ be a subgroup of a group $G$. Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% The set of such right cosets is often written $H \\ G = \{ H g : g \in G \}$.
%
% RepLAB works with right cosets by distinguishing a canonical representative for each coset.
% The set of such canonical representatives is the canonical transversal.

    properties (SetAccess = protected)
        group % (`.NiceFiniteGroup`): Group
        subgroup % (`.NiceFiniteGroup`): Subgroup of `.group`
    end

    methods

        function t = canonicalRepresentative(self, g)
        % Returns the coset representative corresponding to the given element
        %
        % It is thus guaranteed that if ``t = canonicalRepresentative(g)``, then $g t^{-1} = h \in H$, with
        % the decomposition $g = h t$.
        %
        % Moreover ``R.canonicalRepresentative(g) == R.canonicalRepresentative(compose(h, g))`` for any $h \in H$.
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
            C = cellfun(@(h) self.group.compose(h, g), self.subgroup.elements.toCell, 'uniform', 0);
        end

        function C = cosets(self)
        % Returns the set of right cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of cell(1,\*) of elements of `.group`: Set of right cosets
            T = self.transversal;
            C = cellfun(@(t) self.coset(t), T, 'uniform', 0);
        end

    end

end
