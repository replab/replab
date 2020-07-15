classdef NormalCoset < replab.LeftCoset & replab.RightCoset
% Describes a coset in the normal subgroup of a finite group
%
% It is a left and right coset at the same time.

    methods

        function self = NormalCoset(group, canonicalRepresentative)
            self@replab.LeftCoset(group, canonicalRepresentative);
            self@replab.RightCoset(group, canonicalRepresentative);
        end

    end

    methods (Static)

        function c = make(group, element)
            iso = group.niceMorphism;
            chain = iso.image.lexChain;
            % left representative is faster
            permRep = replab.bsgs.Cosets.leftRepresentative(iso.image.lexChain, element);
            c = replab.NormalCoset(group, iso.preimageElement(permRep));
        end

    end

    methods

        function s = cardinality(self)
        % Returns the size of this coset
        %
        % Returns:
        %   vpi: Coset size
            s = cardinality@replab.LeftCoset(self);
        end

        function b = contains(self, el)
        % Returns if this coset contains the given element
        %
        % Args:
        %   el (element of `.type`): Element to check
        %
        % Returns:
        %   logical: True if this coset contains the element
            b = contains@replab.LeftCoset(self, el);
        end

        function E = computeElements(self)
        % Returns an indexed family of the elements of this coset
        %
        % Returns:
        %   `+replab.IndexedFamily`: Elements
            E = computeElements@replab.LeftCoset(self);
        end

    end

end
