classdef NormalCoset < replab.LeftCoset & replab.RightCoset
% Describes a coset in the normal subgroup of a finite group
%
% It is a left and right coset at the same time.
%
% See parent classes `.LeftCoset` and `.RightCoset`

    methods

        % Domain
        function s = sample(self)
            s = sample@replab.LeftCoset(self);
        end

        % FiniteSet

        function C = imap(self, f)
            group1 = self.group.imap(f);
            subgroup1 = self.subgroup.imap(f);
            rep1 = f.imageElement(self.representative);
            C = subgroup1.normalCoset(rep1, 'group', group1);
        end

        function s = nElements(self)
            s = nElements@replab.LeftCoset(self);
        end

        function s = setProduct(self)
            s = setProduct@replab.LeftCoset(self);
        end

    end

end
