classdef DoubleCoset < replab.FiniteSet
% Describes a double coset in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    properties (SetAccess = protected)
        parent % (`.FiniteGroup`): Group containing this double coset
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        H % (`.FiniteGroup`): Subgroup of `.group`
        K % (`.FiniteGroup`): Subgroup of `.group`
        Hchain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
        Kchain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
    end

    methods

        function self = DoubleCoset(H, canonicalRepresentative, K, parent)
            assert(parent.hasSameTypeAs(H));
            assert(parent.hasSameTypeAs(K));
            self.type = parent.type;
            self.parent = parent;
            self.isomorphism = parent.niceMorphism;
            self.H = H;
            self.K = K;
            self.representative = canonicalRepresentative;
            self.Hchain = parent.niceMorphism.imageGroup(H).lexChain;
            self.Kchain = parent.niceMorphism.imageGroup(K).lexChain;
        end

    end

    methods (Static)

        function c = lexCompare(lhs, rhs)
        % Returns the comparison of two permutations by lexicographic order
        %
        % Args:
        %   lhs (permutation): First permutation
        %   rhs (permutation): Second permutation
        %
        % Returns:
        %   integer: 1 if ``lhs > rhs``, 0 if ``lhs == rhs``, -1 if ``lhs < rhs``
            v = lhs - rhs;
            ind = find(v ~= 0, 1)
            c = sign(v(ind));
        end

        function d = make(H, element, K, parent)
            if nargin < 4 || isempty(parent)
                parent = H.closure(K).closure(element);
            end
            canRep = parent.doubleCosets(H, K).cosetRepresentative(element);
            d = replab.DoubleCoset(H, canRep, K, parent);
        end

    end

end
