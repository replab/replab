classdef ConjugacyClass < replab.FiniteSet
% Describes a conjugacy class of a finite group
%
% A conjugacy class containing the representative $r \in G$ is the set $\{g r g^{-1} : g \in G \}$.
%
% The centralizer of $r$ in $G$ is the subgroup $C_{G}(r) = \{ g r g^{-1} == r : g \in G \}$.
%
% Thus, the left cosets $G/C_{G}(r) = \{ g C_{G}(r) : g \in G \}$ are in one to one correspondence with
% the elements of the conjugacy class.
%
% Contrary to the case of cosets, the `.representative` element is not chosen to be lexicographically
% minimal.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group containing this conjugacy class
        representativeCentralizer % (`+replab.FiniteGroup`): Centralizer of `.representative` in `.group`
    end

    methods

        function self = ConjugacyClass(group, arbitraryRepresentative, representativeCentralizer)
            self.type = group.type;
            self.group = group;
            self.representative = arbitraryRepresentative;
            if nargin < 3 || isempty(representativeCentralizer)
                representativeCentralizer = self.group.centralizer(arbitraryRepresentative);
            end
            self.representativeCentralizer = representativeCentralizer;
        end

    end

    methods

        function s = cardinality(self)
            s = self.group.order / self.representativeCentralizer.order;
        end

        function b = contains(self, t)
        % Returns whether the given element is part of this conjugacy class
        %
        % Args:
        %   t (element of `.group`): Group element
        %
        % Returns:
        %   logical: True if ``t`` is a member of this conjugacy classx
            s = self.representative;
            sCentralizer = self.representativeCentralizer;
            % We want to solve ``t == b s b^-1`` with ``s`` the representative
            B = self.group.findLeftConjugations(s, t, sCentralizer);
            b = ~isempty(B);
        end

        function E = computeElements(self)
            T = self.group.leftCosetsOf(self.representativeCentralizer).transversal;
            E = cellfun(@(t) self.group.leftConjugate(t, self.representative), T, 'uniform', 0);
        end

    end

end
