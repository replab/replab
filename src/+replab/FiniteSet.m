classdef FiniteSet < replab.Domain
% Describes a finite set of elements
%
% This is the base class for distinguished subsets of a `.FiniteGroup`, such as cosets or conjugacy classes.
%
% When non-empty, such sets have a distinguished `.representative` element; for permutation groups, this element
% is chosen to be minimal under the lexicographic order.

    properties (SetAccess = protected)
        type % (`.FiniteSet`): Set of all elements of the same type as this set; satisfies ``type.type == type``
        representative % (element of `.type`): Distinguished member of this set; if the set is empty, this value is arbitrary
    end

    methods

        function s = size(self)
        % Returns the size of this set
        %
        % Returns:
        %   vpi: Set cardinality
            error('Abstract');
        end

        function b = contains(self, el)
        % Tests whether this set contains the given element
        %
        % The element must be part of the `.type` set.
        %
        % Args:
        %   el (element of `.type`): Element to test for membership
        %
        % Returns:
        %   logical: True if this set contains ``el`` and false otherwise
            error('Abstract');
        end

        function E = elements(self)
        % Returns an indexed family corresponding to this set
        %
        % Returns:
        %   `.IndexedFamily`: An enumeration of the set elements
            E = self.cached('elements', @() self.computeElements);
        end

        function E = computeElements(self)
        % See `.elements`
            error('Abstract');
        end

    end

    methods % Relations to other sets

        function res = hasSameTypeAs(self, rhs)
        % Returns if this finite set has the same type as the given finite set
        %
        % In particular, it means that the `.contains` method of one set can be called with elements of the other set.
        %x
        % Args:
        %   rhs (`+replab.FiniteSet`): Other finite set
        %
        % Returns:
        %   logical: True if the groups have compatible types
            res = self.type.hasSameTypeAs(rhs.type); % we delegate to the types themselves
        end

    end

end
