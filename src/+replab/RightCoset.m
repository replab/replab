classdef RightCoset < replab.Coset
% Describes a right coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The right coset is $H g = \{ h g : h \in H\}$.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.RightCosetLaws(self);
        end


        % Domain

        function s = sample(self)
            s = self.type.compose(self.subgroup.sample, self.representative);
        end

        % FiniteSet

        function C = imap(self, f)
            group1 = self.group.imap(f);
            subgroup1 = self.subgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = subgroup1.rightCoset(rep1, 'group', group1, 'isCanonical', f.preservesTypeOrder);
        end

    end

end
