classdef FiniteSet < replab.Domain
% Describes a subset of a finite group
%
% Examples of finite subsets include:
%
% - left cosets (`.LeftCoset`),
% - right cosets (`.RightCoset`),
% - normal cosets (`.NormalCoset`),
% - double cosets (`.DoubleCoset`),
% - conjugacy classes (`.ConjugacyClass`),
% - finite groups (`.FiniteGroup`).
%
% If the finite set is not empty (`.nElements` is nonzero), the set has a distinguished `.representative` element.
% For permutation groups, this representative element is the minimal element under lexicographic ordering.

    properties (SetAccess = protected)
        type % (`.FiniteGroupType`): Type of the contained elements
    end

    properties (Access = protected)
        representative_ % (element of `.type`): Minimal member of this set under the `.type` ordering. If the set is empty, value is undefined.
    end

    methods % Elements

        function b = contains(self, el)
        % Tests whether this set contains the given element
        %
        % The element must be part of `.type`, as in the example below: only permutations with the same domain size
        % should be tested using `.contains`.
        %
        % Example:
        %   >>> G = replab.PermutationGroup.of([2 3 4 1]);
        %   >>> G.contains([4 3 2 1])
        %       0
        %
        % Args:
        %   el (element of `.type`): Element to test for membership
        %
        % Returns:
        %   logical: True if this set contains ``el`` and false otherwise
            error('Abstract');
        end

        function E = elements(self)
        % Returns a cell array containing all the elements of this set
        %
        % The sequence is sorted according to the total ordering defined by `.type`.
        %
        % Note: if the number of elements is bigger than `+replab.globals.maxElements`, an error is thrown.
        %
        % Returns:
        %   cell(1,\*) of finite set elements: Elements
            if self.nElements > replab.globals.maxElements
                error('replab:tooBig', 'Number of elements %s too big > %d', num2str(self.nElements), replab.globals.maxElements);
            end
            E = self.elementsSequence.toCell;
        end

        function E = elementsSequence(self)
        % Returns a sequence corresponding to this set
        %
        % The sequence is sorted according to the total ordering defined by `.type`.
        %
        % Returns:
        %   `.Sequence`: An enumeration of the set elements
            error('Abstract');
        end

        function n = nElements(self)
        % Returns the size of this set
        %
        % Returns:
        %   vpi: Set cardinality
            error('Abstract');
        end

        function r = representative(self)
        % Returns the minimal element of this set under the type ordering
        %
        % Returns:
        %   element: Canonical representative of this set
            r = self.representative_;
        end

        function s = setProduct(self)
        % Returns a description of this set as a product of sets
        %
        % Returns:
        %   `.SetProduct`: Generic description of the elements of this set
            error('Abstract');
        end

    end

    methods % Relations to other sets

        function res = hasSameTypeAs(self, rhs)
        % Returns if this finite set has the same type as the given finite set
        %
        % In particular, it means that the `.contains` method of one set can be called with elements of the other set.
        %
        % Args:
        %   rhs (`+replab.FiniteSet`): Other finite set
        %
        % Returns:
        %   logical: True if the groups have compatible types
            res = self.type.isSameTypeAs(rhs.type); % we delegate to the types themselves
        end

    end

    methods % Image under isomorphism

        function res = imap(self, f)
        % Returns the image of this finite set under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism such that this finite set is contained in ``f.source``
        %
        % Returns:
        %   `.FiniteSet`: This finite set mapped under ``f``, expressed as a subset of ``f.image``
            error('Abstract');
        end

    end


    methods % Implementations

        function res = eq(self, rhs)
        % Tests finite sets for equality
        %
        % Example:
        %   >>> G = replab.PermutationGroup.of([2 3 1], [2 1 3]);
        %   >>> H = replab.PermutationGroup.of([2 1 3], [1 3 2]);
        %   >>> G == H
        %       1
            res = replab.finite.equality(self, rhs);
        end

        function res = isequal(self, rhs)
        % Tests finite sets for equality
        %
        % Alternative syntax to `.eq`
        %
        % Example:
        %   >>> G = replab.PermutationGroup.of([2 3 1], [2 1 3]);
        %   >>> H = replab.PermutationGroup.of([2 1 3], [1 3 2]);
        %   >>> isequal(G, H)
        %       1
            res = self == rhs;
        end

        function res = ne(self, rhs)
        % Test finite sets for non-equality
        %
        % Example:
        %   >>> G = replab.PermutationGroup.of([2 3 1], [2 1 3]);
        %   >>> H = replab.PermutationGroup.of([2 1 3]);
        %   >>> G ~= H
        %       1
            res = ~(self == rhs);
        end

        % Domain

        function b = eqv(self, lhs, rhs)
            b = self.type.eqv(lhs, rhs);
        end

        function res = mtimes(self, lhs, rhs)
            res = replab.FiniteSet.multiply(lhs, rhs);
        end

    end

    methods (Static)

        function t = shortType(obj)
            if isa(obj, 'replab.FiniteGroup')
                t = 'G';
            elseif isa(obj, 'replab.ConjugacyClass')
                t = 'C';
            elseif isa(obj, 'replab.DoubleCoset')
                t = 'D';
            elseif isa(obj, 'replab.NormalCoset')
                t = 'N';
            elseif isa(obj, 'replab.LeftCoset')
                t = 'L';
            elseif isa(obj, 'replab.RightCoset')
                t = 'R';
            elseif isa(obj, 'replab.FiniteSet')
                t = 'S';
            else
                t = 'E';
            end
        end

        function res = multiply(lhs, rhs)
            switch [replab.FiniteSet.shortType(lhs) replab.FiniteSet.shortType(rhs)]
              case 'NG'
                res = lhs.group.doubleCoset(lhs.representative, rhs);
              case 'GN'
                res = lhs.doubleCoset(rhs.representative, rhs.group);
              case 'GL'
                res = lhs.doubleCoset(rhs.representative, rhs.group);
              case 'RG'
                res = lhs.group.doubleCoset(lhs.representative, rhs);
              case 'GE'
                res = lhs.rightCoset(rhs);
              case 'EG'
                res = rhs.leftCoset(lhs);
              case 'NE'
                res = lhs.group.rightCoset(lhs.group.type.compose(lhs.representative, rhs));
              case 'EN'
                res = rhs.group.leftCoset(rhs.group.type.compose(lhs, rhs.representative));
              case 'RE'
                res = lhs.group.rightCoset(lhs.group.type.compose(lhs.representative, rhs));
              case 'EL'
                res = rhs.group.leftCoset(rhs.group.type.compose(lhs, rhs.representative));
              otherwise
                error('Invalid product');
                % other cases are invalid
            end
        end

    end

end
