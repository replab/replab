classdef ConjugacyClass < replab.FiniteSet
% Describes a conjugacy class of a finite group
%
% A conjugacy class containing the representative $r \in G$ is the set $\{g r g^{-1} : g \in G \}$.
%
% The centralizer of $r$ in $G$ is the subgroup $C_{G}(r) = \{ g r g^{-1} == r : g \in G \}$.
%
% Thus, the left cosets $G/C_{G}(r) = \{ g C_{G}(r) : g \in G \}$ are in one to one correspondence with
% the elements of the conjugacy class.

    properties (SetAccess = protected) % TODO: Access = protected
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
    end

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group containing this conjugacy class
        representative % (group element): Representative element of the conjugacy class
        representativeCentralizer % (`+replab.FiniteGroup`): Centralizer of `.representative` in `.group`
    end

    methods (Access = protected)

        function self = ConjugacyClass(group, representative, representativeCentralizer)
            self.type = self.group.type;
            self.group = group;
            self.representative = representative;
            self.representativeCentralizer = representativeCentralizer;
        end

    end

    methods

        function s = cardinality(self)
            s = self.group.order / self.representativeCentralizer.order;
        end

        function b = contains(self, el)
        % We want to solve ``el == x c r c^-1`` with ``r`` the representative
        % and ``c`` an element of its centralizer

            b = (self.elements.find(el) ~= 0);
        end

        function E = computeElements(self)
            M = replab.bsgs.Cosets.rightTransversalMatrix(
        end

    end

    methods (Static)

        function c = computeAll(group)
            classes = replab.nfg.conjugacyClassesByOrbits(group);
            c = cellfun(@(cl) replab.ConjugacyClass(group, cl), classes, 'uniform', 0);
        end

    end

end
