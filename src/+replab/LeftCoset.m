classdef LeftCoset < replab.Coset
% Describes a left coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The left coset is $g H = \{ g h : h \in H\}$.

    methods % Implementations

        % Domain

        function s = sample(self)
            s = self.type.compose(self.representative, self.subgroup.sample);
        end

        % FiniteSet

        function C = imap(self, f)
            group1 = self.group.imap(f);
            subgroup1 = self.subgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = subgroup1.leftCoset(rep1, 'group', group1, 'isCanonical', f.preservesTypeOrder);
        end

    end

end
