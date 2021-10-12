classdef FiniteSet < replab.FiniteSet
% Describes a subset of a finite group using a group isomorphism

    properties (SetAccess = protected)
        generic % (`+replab.FiniteSet`): Generic object
        genericIsomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to an object where computations are outsourced
    end

    methods

        function self = FiniteSet(type, generic, genericIsomorphism)
            self.type = type;
            self.generic = generic;
            self.genericIsomorphism = genericIsomorphism;
            self.representative = genericIsomorphism.preimageElement(generic.representative);
        end

    end

    methods

        function b = contains(self, el)
        % Tests whether this set contains the given element
        %
        % The element must be part of ``self.parent.type``.
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

        function s = nElements(self)
        % Returns the size of this set
        %
        % Returns:
        %   vpi: Set cardinality
            error('Abstract');
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
        %x
        % Args:
        %   rhs (`+replab.FiniteSet`): Other finite set
        %
        % Returns:
        %   logical: True if the groups have compatible types
            res = self.type.hasSameTypeAs(rhs.type); % we delegate to the types themselves
        end

    end

    methods % Implementations

        function res = mtimes(lhs, rhs)
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

        function checkNormalizes(lhs, rhs)
            assert(lhs.isNormalizedBy(rhs) || rhs.isNormalizedBy(lhs), 'This product requires a normalization condition');
        end

        function res = multiply(lhs, rhs)
            switch [replab.FiniteSet.shortType(lhs) replab.FiniteSet.shortType(rhs)]
              case 'NG'
                res = replab.DoubleCoset.make(lhs.group, lhs.representative, rhs);
              case 'GN'
                res = replab.DoubleCoset.make(lhs, rhs.representative, rhs.group);
              case 'GL'
                res = replab.DoubleCoset.make(lhs, rhs.representative, rhs.group);
              case 'RG'
                res = replab.DoubleCoset.make(lhs.group, lhs.representative, rhs);
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
