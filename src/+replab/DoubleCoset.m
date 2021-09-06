classdef DoubleCoset < replab.FiniteSet
% Describes a double coset in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    properties (SetAccess = protected)
        parent % (`.FiniteGroup`): Group containing this double coset
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        H % (`.FiniteGroup`): Group
        K % (`.FiniteGroup`): Group
        Hprmgrp % (`.PermutationGroup`): Realization of `.H` as a permutation group
        Kprmgrp % (`.PermutationGroup`): Realization of `.K` as a permutation group
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
            self.Hprmgrp = parent.niceMorphism.imageGroup(H);
            self.Kprmgrp = parent.niceMorphism.imageGroup(K);
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, lhs, rhs)
            b = self.type.eqv(lhs, rhs);
        end

        function l = laws(self)
            l = replab.laws.DoubleCosetLaws(self);
        end

        function s = sample(self)
            s = self.parent.compose(self.H.sample, self.parent.compose(self.representative, self.K.sample));
        end

        % FiniteSet

        function s = size(self)
        % Returns the size of this double coset
        %
        % Returns:
        %   vpi: Coset size

        % From Wikipedia: |H x K| = |H| |K| / |K \intersection x^-1 H x|
            Hconj = self.H.leftConjugateGroup(self.parent.inverse(self.representative));
            inter = self.K.intersection(Hconj);
            s = self.H.order * self.K.order / inter.order;
        end

        function b = contains(self, el)
        % Returns if this double coset contains the given element
        %
        % Args:
        %   el (element of `.type`): Element to check
        %
        % Returns:
        %   logical: True if this coset contains the element
            if ~self.parent.contains(el)
                b = false;
                return
            end
            dc = replab.DoubleCoset.make(self.H, el, self.K, self.parent);
            b = self.parent.eqv(self.representative, dc.representative);
        end

    end

    methods (Access = protected)

        function E = computeElementsSequence(self)
        % Returns an indexed family of the elements of this double coset
        %
        % Returns:
        %   `+replab.Sequence`: Elements
            S = replab.perm.Set(self.Hprmgrp.domainSize);
            Kmat = self.Kprmgrp.chain.allElements;
            g = self.representative';
            S.insert(g(Kmat));
            toCheck = 1;
            while ~isempty(toCheck)
                cur = S.at(toCheck(end))';
                toCheck = toCheck(1:end-1);
                for j = 1:self.Hprmgrp.nGenerators
                    h = self.Hprmgrp.generator(j);
                    newEl = h(cur);
                    if S.find(newEl') == 0
                        newEl = newEl';
                        sz = S.nElements;
                        inds = S.insert(newEl(Kmat));
                        assert(all(inds > sz));
                        toCheck = [toCheck sz+1];
                    end
                end
            end
            S.sort;
            E = replab.seq.FiniteGroupSequence(S.matrix, self.isomorphism);
        end

        function s = computeSetProduct(self)
            s = replab.SetProduct(self.parent, horzcat(Hs.sets, {{self.representative}}, Ks.sets), false);
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
            ind = find(v ~= 0, 1);
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
