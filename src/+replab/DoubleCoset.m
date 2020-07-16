classdef DoubleCoset < replab.FiniteSet
% Describes the set of double cosets in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    properties (SetAccess = protected)
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
        group % (`.FiniteGroup`): Group that includes both `.H` and `.K` as subgroups
        H % (`.FiniteGroup`): Subgroup of `.group`
        K % (`.FiniteGroup`): Subgroup of `.group`
        Hchain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
        Kchain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
        leftCosets % (`+replab.LeftCosets`): Left cosets G / K
    end

    methods

        function self = DoubleCosets(group, H, K)
            assert(group.hasSameTypeAs(H));
            assert(group.hasSameTypeAs(K));
            self.group = group;
            self.H = H;
            self.K = K;
            self.isomorphism = group.niceMorphism;
            self.groupChain = self.isomorphism.imageGroup(group).lexChain;
            self.Hchain = self.isomorphism.imageGroup(H).lexChain;
            self.Kchain = self.isomorphism.imageGroup(K).lexChain;
            self.leftCosets = group / K;
        end

    end

end
