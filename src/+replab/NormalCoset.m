classdef NormalCoset < replab.LeftCoset & replab.RightCoset
% Describes a coset in the normal subgroup of a finite group
%
% It is a left and right coset at the same time.
%
% See parent classes `.LeftCoset` and `.RightCoset`

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.NormalCosetLaws(self);
        end

        % Domain

        function s = sample(self)
            s = sample@replab.LeftCoset(self);
        end

        % FiniteSet

        function C = imap(self, f)
            group1 = self.group.imap(f);
            subgroup1 = self.subgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = subgroup1.normalCoset(rep1, 'group', group1, 'isCanonical', f.preservesTypeOrder);
        end

    end

end
