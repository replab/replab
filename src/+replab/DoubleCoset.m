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
            Hprmgrp = parent.niceMorphism.imageGroup(H);
            Kprmgrp = parent.niceMorphism.imageGroup(K);
            S = replab.perm.Set(Hprmgrp.domainSize);
            elP = parent.niceMorphism.imageElement(element);
            repP = replab.bsgs.Cosets.leftRepresentative(Kprmgrp.lexChain, elP);
            minP = repP;
            S.insert(repP');
            toCheck = 1;
            % compute the orbit of the left cosets ``element K``
            while ~isempty(toCheck)
                i = toCheck(end);
                toCheck = toCheck(1:end-1);
                rep = S.at(i)';
                for j = 1:Hprmgrp.nGenerators
                    h = Hprmgrp.generator(j);
                    elP = h(rep);
                    repP = replab.bsgs.Cosets.leftRepresentative(Kprmgrp.lexChain, elP);
                    ind = S.find(repP');
                    if ind == 0
                        if self.lexCompare(repP, minP) < 0
                            minP = repP;
                        end
                        ind = S.insert(repP');
                        toCheck = [toCheck ind];
                    end
                end
            end
            d = replab.DoubleCoset(H, parent.niceMorphism.preimageElement(minP), K, parent);
        end

    end

end
