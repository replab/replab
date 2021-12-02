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

        function l = laws(self)
            l = replab.laws.DoubleCosetLaws(self);
        end

        function s = sample(self)
            s = self.type.compose(self.leftSubgroup.sample, self.type.compose(self.representative, self.rightSubgroup.sample));
        end

        % FiniteSet

        function C = imap(self, f)
            group1 = self.group.imap(f);
            leftSubgroup1 = self.leftSubgroup.imap(f);
            rightSubgroup1 = self.rightSubgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = leftSubgroup1.doubleCoset(rep1, rightSubgroup1, 'group', group1, 'isCanonical', f.preservesTypeOrder);
        end

    end

end
