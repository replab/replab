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

        function s = nElements(self)
            s = nElements@replab.LeftCoset(self);
        end

        function s = setProduct(self)
            s = setProduct@replab.LeftCoset(self);
        end

    end

end
