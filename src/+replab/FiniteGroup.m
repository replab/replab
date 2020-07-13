classdef FiniteGroup < replab.CompactGroup & replab.FiniteSet
% Describes a group with a finite number of elements
%
% Each finite group has a type, that describes the most general group embedding its elements. For example, permutations of domain size ``n`` are embedded in the symmetric group of degree ``n``.

    properties (SetAccess = protected)
        generators % (cell(1,\*) of `.type` elements): Group generators
    end

    methods % Implementations

        function res = ne(self, rhs)
            res = ~(self == rhs);
        end

        function res = eq(self, rhs)
            res = self.hasSameTypeAs(rhs) && self.isSubgroupOf(rhs) && rhs.isSubgroupOf(self);
        end

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'generators';
            names{1, end+1} = 'representative';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nGenerators
                names{1, end+1} = sprintf('generator(%d)', i);
                values{1, end+1} = self.generator(i);
            end
        end

    end

    methods % Group properties

        function s = cardinality(self)
            s = self.order;
        end

        function o = order(self)
        % Returns the group order
        %
        % Returns:
        %   vpi: The group order
            o = self.cached('order', @() self.computeOrder);
        end

        function o = computeOrder(self)
        % See `.order`
            error('Abstract');
        end

        function D = decomposition(self)
        % Returns a decomposition of this group as a product of sets
        %
        % Returns:
        %   `.FiniteGroupDecomposition`: The group decomposition
            D = self.cached('decomposition', @() self.computeDecomposition);
        end

        function D = computeDecomposition(self)
        % See `.decomposition`
            error('Abstract');
        end

        function e = exponent(self)
        % Returns the group exponent
        %
        % The group exponent is the smallest integer ``e`` such that ``g^e == identity`` for all ``g`` in ``G``.
        %
        % Returns:
        %   vpi: The group exponent
            e = self.cached('exponent', @() self.computeExponent);
        end

        function e = computeExponent(self)
            eo = cellfun(@(c) self.elementOrder(c.representative), self.conjugacyClasses);
            eo = unique(eo);
            e = eo(1);
            for i = 2:length(eo)
                e = lcm(e, eo(i));
            end
        end

        function c = conjugacyClasses(self)
        % Returns the conjugacy classes of this group
        %
        % Returns:
        %   cell(1, \*) of `+replab.ConjugacyClass`: Array of conjugacy classes
            c = self.cached('conjugacyClasses', @() self.computeConjugacyClasses);
        end

        function c = computeConjugacyClasses(self)
            error('Abstract');
        end

        function R = recognize(self)
        % Attempts to recognize this group in the standard atlas
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case the group is identified; or ``[]`` if unrecognized.
            R = self.cached('recognize', @() self.computeRecognize);
        end

        function R = computeRecognize(self)
            A = replab.atlas.Standard;
            R = A.recognize(self);
        end

        function b = isTrivial(self)
        % Tests whether this group is trivial
        %
        % Returns:
        %   logical: True if this group is trivial (i.e. has only one element)
            b = self.nGenerators == 0;
        end

        function res = isCommutative(self)
        % Returns whether this group is commutative
            for i = 1:self.nGenerators
                gi = self.generator(i);
                for j = 1:i-1
                    gj = self.generator(j);
                    if ~self.eqv(self.compose(gi, gj), self.compose(gj, gi))
                        res = false;
                        return
                    end
                end
            end
            res = true;
            return
        end

        function res = isCyclic(self)
        % Returns whether this group is a cyclic group
            res = self.cached('isCyclic', @() self.computeIsCyclic);
        end

        function res = computeIsCyclic(self)
            error('Abstract');
        end

        function res = isSimple(self)
        % Returns whether this group is simple
            res = self.cached('isSimple', @() self.computeIsSimple);
        end

        function res = computeIsSimple(self)
            error('Abstract');
        end

    end

    methods % Group elements

        function n = nGenerators(self)
        % Returns the number of group generators
        %
        % Returns:
        %   integer: Number of group generators
            n = length(self.generators);
        end

        function p = generator(self, i)
        % Returns the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: i-th group generator
            p = self.generators{i};
        end

        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: Inverse of the i-th group generator
            p = self.inverse(self.generators{i});
        end

        function b = contains(self, g)
        % Tests whether this group contains the given element
        %
        % Args:
        %   g (element of `.type`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains ``g`` and false otherwise
            error('Abstract');
        end

        function o = elementOrder(self, g)
        % Returns the order of a group element
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   vpi: The order of ``g``, i.e. the smallest ``o`` such that ``g^o == identity``
            error('Abstract');
        end

    end

    methods % Construction of groups

        function conj = leftConjugateGroup(self, by)
        % Returns the left conjugate of the current group by the given element
        %
        % ``res = self.leftConjugateGroup(by)``
        %
        % In particular, it ensures that
        % ``res.generator(i) = self.type.leftConjugate(by, self.generator(i))``
        %
        % Args:
        %   by (element of `.type`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.FiniteGroup`: The conjugated group
            newGenerators = cellfun(@(g) self.type.leftConjugate(by, g), self.generators, 'uniform', 0);
            conj = self.type.subgroupWithGenerators(newGenerators, self.order);
        end

        function res = closure(self, obj)
        % Computes the group generated by the elements of this group and another object
        %
        % Example:
        %    >>> G = replab.CyclicGroup(3);
        %    >>> S3 = replab.S(3);
        %    >>> H = S3.subgroup({[2 1 3]});
        %    >>> C = G.closure(H);
        %    >>> C == S3
        %        1
        %
        % Args:
        %   rhs (`.FiniteGroup` or group element): Element(s) to add
        %
        % Returns:
        %   `.FiniteGroup`: The closure
            error('Abstract');
        end

        function res = normalClosure(self, obj)
        % Computes the normal closure of an object in the closure of this group and this object
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> S4 = S5.subgroup({[2 3 4 1 5] [2 1 3 4 5]});
        %   >>> nc = S4.normalClosure(S5.subgroup({[1 2 4 5 3]}));
        %   >>> nc == S5.subgroup({[1 2 4 5 3] [5 2 3 1 4] [2 5 3 4 1]})
        %       1
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> nc = S4.normalClosure(S4.subgroup({[2 3 1 4]}));
        %   >>> nc == S4.subgroup({[2 3 1 4] [1 3 4 2]});
        %       1
        %
        % Args:
        %   obj (`+replab.FiniteGroup` or group element): Element to compute the normal closure with
        %
        % Returns:
        %   `+replab.FiniteGroup`: The normal closure of ``obj`` in this group
            error('Abstract');
        end

    end

    methods % Subgroups

        function sub = subgroupWithGenerators(self, generators, order)
        % Constructs a subgroup of the current group from generators
        %
        % Guarantees that the constructed subgroup preserves the order of the generators. However, the
        % given generators cannot contain the identity element.
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> G = S5.subgroupWithGenerators({[2 3 4 5 1]});
        %   >>> G.generators
        %       TODO
        %
        % Args:
        %   generators (cell(1,\*) of elements of this group): List of generators, does not contain the identity
        %   order (vpi or ``[]``, optional): Subgroup order
        %
        % Returns:
        %   `+replab.FiniteGroup`: The subgroup generated by `.generators`
            error('Abstract');
        end

        function sub = subgroup(self, elements, order)
        % Constructs a subgroup of the current group from elements
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> G = S5.subgroupWithGenerators({[1 2 3 4 5]});
        %   >>> G.generators
        %       TODO is empty
        %
        % Args:
        %   elements (cell(1,\*) of elements of this group): List of elements
        %   order (vpi, optional): Subgroup order
        %
        % Returns:
        %   `+replab.FiniteGroup`: The subgroup generated by `.generators`
            if nargin < 3
                order = [];
            end
            mask = cellfun(@(e) self.isIdentity(e), elements);
            sub = self.subgroupWithGenerators(elements(~mask), order);
        end

        function sub = trivialSubgroup(self)
        % Returns the trivial subgroup of this group
        %
        % Returns:
        %   `+replab.FiniteGroup`: The trivial subgroup
            sub = self.subgroup({}, vpi(1));
        end

        function sub = derivedSubgroup(self)
        % Computes the derived subgroup of this group
        %
        % Example:
        %   >>> S4 = replab.SymmetricGroup(4);
        %   >>> S4.order
        %       24
        %   >>> D = S4.derivedSubgroup;
        %   >>> D == replab.AlternatingGroup(4)
        %       1
        % Returns:
        %   `+replab.FiniteGroup`: The derived subgroup
            sub = self.cached('derivedSubgroup', @() self.computeDerivedSubgroup);
        end

        function sub = computeDerivedSubgroup(self)
        % See `.derivedSubgroup`
            error('Abstract');
        end

        function sub = center(self)
        % Returns the center of this group
        %
        % Returns:
        %   `.FiniteGroup`: The center
            sub = self.cached('center', @() self.centralizer(self));
        end

        function sub = centralizer(self, arg)
        % Returns the centralizer of the given object in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizer([2 3 1 4]);
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizer(G.subgroup({[2 3 1 4]}));
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Args:
        %   rhs (`.FiniteGroup` or group element): Element to compute the centralizer of
        %
        % Returns:
        %   `.FiniteGroup`: The centralizer
            error('Abstract');
        end

        function sub = intersection(self, G)
        % Computes the intersection of two groups
        %
        % Example:
        %   >>> D12 = replab.S(6).subgroup({[1 6 5 4 3 2], [2 1 6 5 4 3]});
        %   >>> G = replab.S(6).subgroup({[2 1 3 4 5 6], [2 3 4 5 1 6]});
        %   >>> I = D12.intersection(G);
        %   >>> I == replab.S(6).subgroup({[5 4 3 2 1 6]})
        %       1
        %
        % Args:
        %   G (`+replab.FiniteGroup`): Group with the same type as this group
        %
        % Returns:
        %   `+replab.FiniteGroup`: Subgroup representing the intersection
            error('Abstract');
        end

    end

    methods % Cosets

        function c = normalCoset(self, normalSubgroup, element)
        % Returns a normal coset
        %
        % Returns the set ``element * normalSubgroup == normalSubgroup * element``.
        %
        % Args:
        %   normalSubgroup (`+replab.FiniteGroup`): Normal subgroup of this group
        %   element (group element): Group element
        %
        % Returns:
        %   `+replab.NormalCoset`: The constructed normal coset
            c = self.normalCosetsOf(subgroup).coset(element);
        end

        function c = normalCosetsOf(self, subgroup)
        % Returns the set of normal cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Normal subgroup of this group
            assert(subgroup.isNormalSubgroupOf(self), 'The given subgroup must be normal in parent group');
            c = replab.NormalCosets(self, subgroup);
        end

        function c = rightCoset(self, subgroup, element)
        % Returns a right coset
        %
        % Returns the set ``subgroup * element``.
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %   element (group element): Group element
        %
        % Returns:
        %   `+replab.RightCoset`: The constructed right coset
            c = self.rightCosetsOf(subgroup).coset(element);
        end

        function c = rightCosetsOf(self, subgroup)
        % Returns the set of right cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %
        % Returns:
        %   `+replab.RightCosets`: Right cosets
            c = replab.RightCosets(self, subgroup);
        end

        function c = mldivide(self, supergroup)
        % Shorthand for `.rightCosetsOf`
            c = supergroup.rightCosetsOf(self);
        end

        function c = leftCoset(self, subgroup, element)
        % Returns a left coset
        %
        % Returns the set ``element * subgroup``.
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %   element (group element): Group element
        %
        % Returns:
        %   `+replab.LeftCoset`: The constructed right coset
            c = self.leftCosetsOf(subgroup).coset(element);
        end

        function c = leftCosetsOf(self, subgroup)
        % Returns the set of left cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %
        % Returns:
        %   `+replab.LeftCosets`: Left cosets
            c = replab.LeftCosets(self, subgroup);
        end

        function c = mrdivide(self, subgroup)
        % Shorthand for `.leftCosetsOf`
            c = self.leftCosetsOf(subgroup);
        end

        function B = findLeftConjugations(self, s, t, sCentralizer, tCentralizer)
        % Returns the set of all elements that left conjugates an element to another element
        %
        % Let ``s`` and ``t`` be two elements of this group. This returns the set of all elements
        % ``b`` such that ``t = b s b^-1`` or ``t = leftConjugate(b, s)``.
        %
        % When no such ``b`` exists, this returns ``[]``.
        %
        % Args:
        %   s (group element): Source element
        %   t (group element): Target element
        %   sCentralizer (`+replab.FiniteGroup` or ``[]``, optional): Centralizer of ``s`` in this group
        %   tCentralizer (`+replab.FiniteGroup` or ``[]``, optional): Centralizer of ``t`` in this group
        %
        % Returns:
        %   `+replab.LeftCoset` or ``[]``: Set of all elements of this group that left conjugates ``s`` to ``t`` if it exists
            error('Abstract');
        end

    end

    methods % Relations to other groups

        function res = isSubgroupOf(self, rhs)
        % Returns whether this group is a subgroup of another group
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> G = S3.subgroup({[2 1 3]});
        %   >>> G.isSubgroupOf(S3)
        %       1
        %
        % Args:
        %   rhs (`+replab.FiniteGroup`): Other group with the same type as this one
        %
        % Returns:
        %   logical: True if this group is a subgroup of ``rhs``
            res = all(cellfun(@(g) rhs.contains(g), self.generators));
        end

        function res = isNormalSubgroupOf(self, rhs)
        % Returns whether this group is a normal subgroup of another group
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> A3 = replab.AlternatingGroup(3);
        %   >>> A3.isNormalSubgroupOf(S3)
        %       1
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> G = S3.subgroup({[2 1 3]});
        %   >>> G.isNormalSubgroupOf(S3)
        %       0
        %
        % Args:
        %   rhs (`+replab.FiniteGroup`): Other group with the same type as this one
        %
        % Returns:
        %   logical: True if this group is a normal subgroup of ``rhs``
            if ~self.isSubgroupOf(rhs)
                res = false;
                return
            end
            for i = 1:self.nGenerators
                subi = self.generator(i);
                for j = 1:rhs.nGenerators
                    rhsj = rhs.generator(j);
                    if ~self.contains(self.leftConjugate(rhsj, subi))
                        res = false;
                        return
                    end
                end
            end
            res = true;
        end

    end

    methods % Morphisms

% $$$         % TODO
% $$$         function f = leftConjugateMorphism(self, by)
% $$$         % Returns the morphism that corresponds to left conjugation by an element
% $$$         %
% $$$         % Args:
% $$$         %   by (element of `parent`): Element to conjugate the group with
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Morphism`: Conjugation morphism
% $$$             generatorImages = cellfun(@(g) self.parent.leftConjugate(by, g), self.generators, 'uniform', 0);
% $$$             target = self.parent.subgroup(newGenerators, self.order);
% $$$             f = self.morphismByImages(self, target, generatorImages);
% $$$         end

        function m = morphismByImages(self, target, generatorImages)
        % Constructs a morphism to a group using images of generators
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> m = S4.morphismByImages(replab.S(3), {[1 3 2] [3 2 1]});
        %   >>> replab.FiniteMorphismLaws(m).checkSilent
        %       1
        %
        % Args:
        %   target (`.FiniteGroup`): Target of the morphism, the morphism image is a subgroup of this
        %   generatorImages (cell(1, \*) of target elements): Images of this group generators
        %
        % Returns:
        %   `.FiniteMorphism`: The constructed morphism
            error('Abstract');
        end

    end

% $$$
% $$$
% $$$     methods % Representations
% $$$
% $$$         function rep = regularRep(self)
% $$$         % Returns the left regular representation of this group
% $$$         %
% $$$         % Returns:
% $$$         %   `.Rep`: The left regular representation as a real permutation representation
% $$$             error('Abstract');
% $$$         end
% $$$
% $$$         function rho = repByImages(self, field, dimension, images)
% $$$         % Constructs a finite dimensional representation of this group from generator images
% $$$         %
% $$$         % Args:
% $$$         %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
% $$$         %   dimension (integer): Representation dimension
% $$$         %   images (cell(1,\*) of double(\*,\*), may be sparse): Images of the group generators
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The constructed group representation
% $$$             error('Abstract');
% $$$         end
% $$$
% $$$         function rho = permutationRep(self, dimension, permutations)
% $$$         % Constructs a permutation representation of this group
% $$$         %
% $$$         % The returned representation is real. Use `+replab.Rep.complexification` to obtain a complex representation.
% $$$         %
% $$$         % Args:
% $$$         %   dimension (integer): Dimension of the representation
% $$$         %   permutations (cell(1,\*) of permutations): Images of the generators as permutations of size ``dimension``
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The constructed group representation
% $$$             f = @(g) replab.Permutation.toSparseMatrix(g);
% $$$             images = cellfun(f, permutations, 'uniform', 0);
% $$$             inverseImages = cellfun(@(i) i', images, 'uniform', 0);
% $$$             rho = self.repByImages('R', dimension, images, inverseImages);
% $$$         end
% $$$
% $$$         function rho = signedPermutationRep(self, dimension, signedPermutations)
% $$$         % Returns a real signed permutation representation of this group
% $$$         %
% $$$         % The returned representation is real. Use ``rep.complexification`` to obtain a complex representation.
% $$$         %
% $$$         % Args:
% $$$         %   dimension: Dimension of the representation
% $$$         %   signedPermutations (cell(1,\*) of signed permutations): Images of the generators as signed permutations of size ``dimension``
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The constructed group representation
% $$$             f = @(g) replab.SignedPermutation.toSparseMatrix(g);
% $$$             images = cellfun(f, signedPermutations, 'uniform', 0);
% $$$             rho = self.repByImages('R', dimension, images);
% $$$         end
% $$$
% $$$     end

% TODO: coset stuff
end
