classdef DoubleCoset < replab.FiniteSet
% Describes a double coset in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group containing this double coset
        leftSubgroup % (`.FiniteGroup`): Group
        rightSubgroup % (`.FiniteGroup`): Group
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
            s = self.type.compose(self.leftSubgroup.sample, self.type.compose(self.representative, self.rightSubgroup.sample));
        end

        % FiniteSet

        function b = contains(self, el)
            if ~self.group.contains(el)
                b = false;
                return
            end
            dc = self.leftSubgroup.doubleCoset(el, self.rightSubgroup, 'group', self.group);
            b = self.type.eqv(self.representative, dc.representative);
        end

        function C = imap(self, f)
            group1 = self.group.imap(f);
            leftSubgroup1 = self.leftSubgroup.imap(f);
            rightSubgroup1 = self.rightSubgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = leftSubgroup1.doubleCoset(rep1, rightSubgroup1, 'group', group1, 'isCanonical', f.preservesTypeOrder);
        end

        function s = nElements(self)

        % From Wikipedia: |left x right| = |left| |right| / |right \intersection x^-1 left x|
            leftConj = self.leftSubgroup.leftConjugateGroup(self.group.inverse(self.representative));
            inter = self.rightSubgroup.intersection(leftConj);
            s = self.leftSubgroup.order * self.rightSubgroup.order / inter.order;
        end

        function s = setProduct(self)
            s = replab.SetProduct(self.type, horzcat(self.leftSubgroup.setProduct.sets, {{self.representative}}, self.rightSubgroup.setProduct.sets), false);
        end

    end

end
