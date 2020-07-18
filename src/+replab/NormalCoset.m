classdef NormalCoset < replab.LeftCoset & replab.RightCoset
% Describes a coset in the normal subgroup of a finite group
%
% It is a left and right coset at the same time.

    methods

        function self = NormalCoset(group, canonicalRepresentative, parent)
            self@replab.LeftCoset(group, canonicalRepresentative, parent);
            self@replab.RightCoset(group, canonicalRepresentative, parent);
        end

    end

    methods (Static)

        function c = make(group, element, parent)
            if nargin < 3 || isempty(parent)
                parent = group.closure(element);
            end
            chain = parent.niceMorphism.imageGroup(group).lexChain;
            permRep = replab.bsgs.Cosets.leftRepresentative(chain, parent.niceMorphism.imageElement(element));
            % left representative is faster
            c = replab.NormalCoset(group, parent.niceMorphism.preimageElement(permRep), parent);
        end

    end

    methods

        function s = sample(self)
            s = sample@replab.LeftCoset(self);
        end

        function s = size(self)
            s = size@replab.LeftCoset(self);
        end

        function b = contains(self, el)
            b = contains@replab.LeftCoset(self, el);
        end

        function E = computeElements(self)
            E = computeElements@replab.LeftCoset(self);
        end

    end

end
