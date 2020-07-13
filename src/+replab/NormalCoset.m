classdef NormalCoset < replab.LeftCoset & replab.RightCoset
% Describes a coset in the normal subgroup of a finite group
%
% It is a left and right coset at the same time.

    methods (Access = protected)

        function self = NormalCoset(group, subgroup, representative)
            self@replab.LeftCoset(group, subgroup, representative);
            self@replab.RightCoset(group, subgroup, representative);
        end

    end

    methods

        function s = size(self)
        % Returns the size of this coset
        %
        % Returns:
        %   vpi: Coset size
            s = size@replab.LeftCoset(self);
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

    methods (Static)

        function L = make(group, subgroup, element)
            assert(isa(group, 'replab.PermutationGroup')); % TODO
            L = replab.NormalCoset(group, subgroup, replab.bsgs.Cosets.leftRepresentative(subgroup.lexChain, element));
        end

        function L = fromRepresentative(group, subgroup, representative)
            L = replab.NormalCoset(group, subgroup, representative);
        end

    end

end
