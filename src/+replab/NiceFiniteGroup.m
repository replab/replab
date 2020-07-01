classdef NiceFiniteGroup < replab.FiniteGroup
% A nice finite group is a finite group equipped with an injective homomorphism into a permutation group
%
% The class that subclasses `.NiceFiniteGroup` implements a method `.niceMonomorphismImage` that returns a
% permutation row vector corresponding to a group element.
%
% In turn, the `.NiceFiniteGroup` infrastructure will use that method to build a BSGS chain to describe
% the structure of the finite group; this chain also provides a way to compute the preimage of a permutation.
%
% Thus, an isomorphism is established between the present `.NiceFiniteGroup` and a permutation group; as
% permutation groups can be handled by efficient BSGS algorithms, the requested computations can be
% translated back and forth between this group and a permutation group.
%
% In particular, the decomposition of the finite group in a product of sets (`.FiniteGroupDecomposition`),
% the enumeration of elements using a `.IndexedFamily`, the construction of subgroups is all handled
% by permutation group algorithms.
%
% Each nice finite group has a parent object, that describes the most general group embedding
% elements of a particular type. For example, permutations of domain size ``n`` are embedded in the
% symmetric group of degree ``n`` (for such groups, their nice monomorphism is the identity).
%
% The `.niceMonomorphismImage` method of this group should be valid for all elements of the parent as well;
% then the `.contains` method will be valid for all elements of the parent group.

    properties (SetAccess = protected)
        parent % `+replab.NiceFiniteGroup`: Parent nice finite group
    end


    properties (Access = protected)
        niceGroup_ % `+replab.+PermutationGroup`: Image of this group through the nice monomorphism
        niceInverseMonomorphism_ % `+replab.mrp.PermMorphism`: Inverse of the nice monomorphism
    end

    methods

        %% Abstract

        function res = hasSameParentAs(self, rhs)
        % Returns whether this group has a parent compatible with another group
        %
        % This method is used to test for group equality
        %
        % Args:
        %   rhs (`+replab.NiceFiniteGroup`): Other nice finite group
        %
        % Returns:
        %   logical: True if the groups have compatible parents
            error('Abstract');
        end

        function p = niceMonomorphismImage(self, g)
        % Returns a permutation representation of the given group element
        %
        % A nice monomorphism is the GAP System terminology for injective
        % homomorphism into a permutation group.
        %
        % Args:
        %   g (element): Group element to represent as a permutation
        %
        % Returns:
        %   permutation: Permutation representation of ``g``
            error('Abstract');
        end

        function P = niceMonomorphismGroupImage(self, G)
        % Returns the image of a subgroup of this nice finite group through the nice monomorphism
        %
        % Args:
        %   G (`+replab.NiceFiniteGroup`): Preimage group, must be a subgroup of this group
        %
        % Returns:
        %   `+replab.PermutationGroup`: Image
            n = length(self.niceMonomorphismImage(self.identity));
            % possible optimization: reuse the stabilizer chain in the inverse morphism of ``G``
            generators = cellfun(@(x) self.niceMonomorphismImage(x), G.generators, 'uniform', 0);
            chain = replab.bsgs.Chain(n, generators);
            P = replab.PermutationGroup(n, gens, chain.order, [], chain);
        end

        function G = niceMonomorphismGroupPreimage(self, P)
        % Returns the preimage of a permutation subgroup of the nice group of this group
        %
        % Args:
        %   P (`+replab.PermutationGroup`): Image permutation group, must be a subgroup of `.niceGroup`
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: Preimage which is a subgroup of this group
            generators = cellfun(@(x) self.niceMonomorphismPreimage(x), P.generators, 'uniform', 0);
            G = self.subgroup(generators);
        end

    end

    methods (Access = protected)

        %% FiniteGroup methods

        function order = computeOrder(self)
            order = self.niceGroup.chain.order;
        end

        function m = computeNiceInverseMonomorphism(self)
            m = replab.mrp.PermMorphism.byImages(self.niceGroup, self, self.generators);
        end

        function g = computeNiceGroup(self)
            imgId = self.niceMonomorphismImage(self.identity);
            n = length(imgId);
            g = replab.PermutationGroup(n, self.niceGenerators, self.order_);
        end

        function E = computeElements(self)
            basis = replab.util.MixedRadix(self.niceGroup.chain.orbitSizes, true, true);
            atFun = @(ind) self.niceMonomorphismPreimage(self.niceGroup.chain.elementFromIndices(basis.sub2ind(ind)));
            findFun = @(el) basis.ind2sub(self.niceGroup.chain.indicesFromElement(self.niceMonomorphismImage(el)));
            E = replab.IndexedFamily.lambda(self.order, atFun, findFun);
        end

        function dec = computeDecomposition(self)
            dec = replab.FiniteGroupDecomposition(self, self.niceInverseMonomorphism.chain.imagesDecomposition);
        end

    end

    methods

        function res = eq(self, rhs)
            res = self.hasSameParentAs(rhs) && self.isSubgroupOf(rhs) && rhs.isSubgroupOf(self);
        end

        function res = isSubgroupOf(self, rhs)
        % Returns whether this group is a subgroup of another group
        %
        % Args:
        %   rhs (`+replab.NiceFiniteGroup`): Another nice finite group with the same parent
        %
        % Returns:
        %   logical: True if this group is a subgroup of ``res``
            res = all(cellfun(@(g) rhs.contains(g), self.generators));
        end

        function sub = subgroup(self, generators, order)
        % Constructs a subgroup of the current group from generators
        %
        % Args:
        %   generators (row cell array of elements of this group): List of generators
        %   order (vpi, optional): Subgroup order
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: The subgroup generated by `.generators`
            if nargin < 3
                order = [];
            end
            sub = replab.NiceFiniteSubgroup(self.parent, generators, order);
        end

        function conj = leftConjugateGroup(self, by)
        % Returns the left conjugate of the current group by the given element
        %
        % ``res = self.leftConjugateGroup(by)``
        %
        % In particular, it ensures that
        % ``res.generator(i) = self.parent.leftConjugate(by, self.generator(i))``
        %
        % Args:
        %   by (element of `parent`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: The conjugated group
            newGenerators = cellfun(@(g) self.leftConjugate(by, g), self.generators, 'uniform', 0);
            conj = self.subgroup(newGenerators, self.order);
        end

        function sub = normalClosure(self, H)
        % Computes the normal closure of a group H in the closure of this group and H
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
        %   rhs (`+replab.NiceFiniteGroup`): A permutation group acting on the same domain as this group
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: The normal closure of ``rhs`` in this group
            assert(~isa(self, 'replab.PermutationGroup'));
            assert(self.hasSameParentAs(H));
            sub = self.parent.niceMonomorphismGroupPreimage(self.niceGroup.normalClosure(H.niceGroup));
        end

        function sub = derivedSubgroup(self)
        % Computes the derived subgroup of this group
        %
        % Example:
        %   >>> S4 = replab.Permutations(4);
        %   >>> S4.order
        %       24
        %   >>> D = S4.derivedSubgroup;
        %   >>> D.order % the derived subgroup of S4 is the alternating group A4
        %       12
        % Returns:
        %   `+replab.NiceFiniteGroup`: The derived subgroup
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            sub = self.niceMonomorphismGroupPreimage(self.niceGroup.derivedSubgroup);
        end

        function grp = trivialSubgroup(self)
        % Returns the trivial subgroup of this group
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: The trivial subgroup
            grp = self.subgroup({}, vpi(1));
        end

        function o = elementOrder(self, g)
        % Returns the order of a group element
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   integer: The order of ``g``, i.e. the smallest ``o`` such that ``g^o == identity``
            o = self.niceGroup.elementOrder(self.niceMonomorphismImage(g));
        end

        %% Methods enabled by the BSGS algorithms

        function c = niceGroup(self)
        % Returns the image of this group as a permutation group through the nice monomorphismx
        %
        % Returns:
        %   `+replab.PermutationGroup`: Permutation group
            if isempty(self.niceGroup_)
                self.niceGroup_ = self.computeNiceGroup;
            end
            c = self.niceGroup_;
        end

        function m = niceInverseMonomorphism(self)
        % Returns the monomorphism from the permutation representation to the original group
            if isempty(self.niceInverseMonomorphism_)
                self.niceInverseMonomorphism_ = self.computeNiceInverseMonomorphism;
            end
            m = self.niceInverseMonomorphism_;
        end

        function ng = niceGenerators(self)
        % Returns the image of the group generators under the nice monomorphism
            ng = cellfun(@(g) self.niceMonomorphismImage(g), self.generators, 'uniform', 0);
        end

        function m = niceMonomorphism(self)
            m = replab.Morphism.lambda(self, self.niceGroup, @(g) self.niceMonomorphismImage(g));
        end

        function g = niceMonomorphismPreimage(self, p)
        % Returns the group element corresponding to a permutation
        %
        % See also `.niceMonomorphismImage`
        %
        % Args:
        %   p (permutation): Permutation representation
        %
        % Returns:
        %   g (element): Group element corresponding to the permutation
            g = self.niceInverseMonomorphism.image(p);
        end

        %% CompactGroup methods

        function g = sampleUniformly(self)
            [~, g] = self.niceInverseMonomorphism.chain.sampleUniformly;
        end


        function b = contains(self, g)
        % Tests whether this group contains the given parent group element
        %
        % Abstract in `+replab.NiceFiniteSubgroup`
        %
        % Args:
        %   g (element of `parent`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains ``g`` and false otherwise
            b = self.niceGroup.chain.contains(self.niceMonomorphismImage(g));
        end

        %% Representation construction

        function rho = repByImages(self, field, dimension, images, inverseImages)
        % Constructs a finite dimensional representation of this group from generator images
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   images (cell(1,\*) of double(\*,\*), may be sparse): Images of the group generators
        %   inverseImages (cell(1,\*) of double(\*,\*), may be sparse): Inverse images of the group generators
        % Returns:
        %   `+replab.Rep`: The constructed group representation
            rho = replab.RepByImages(self, field, dimension, images, inverseImages);
        end

        function rho = permutationRep(self, dimension, permutations)
        % Constructs a permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification`` to obtain a complex representation.
        %
        % Args:
        %   dimension (integer): Dimension of the representation
        %   permutations (cell(1,\*) of permutations): Images of the generators as permutations of size ``dimension``
        %
        % Returns:
        %   `+replab.Rep`: The constructed group representation
            S = replab.Permutations(dimension);
            f = @(g) S.toSparseMatrix(g);
            images = cellfun(f, permutations, 'uniform', 0);
            inverseImages = cellfun(@(i) i', images, 'uniform', 0);
            rho = self.repByImages('R', dimension, images, inverseImages);
        end

        function rho = signedPermutationRep(self, dimension, signedPermutations)
        % Returns a real signed permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification`` to obtain a complex representation.
        %
        % Args:
        %   dimension: Dimension of the representation
        %   signedPermutations (cell(1,\*) of signed permutations): Images of the generators as signed permutations of size ``dimension``
        %
        % Returns:
        %   `+replab.Rep`: The constructed group representation
            f = @(g) replab.SignedPermutations.toSparseMatrix(g);
            images = cellfun(f, signedPermutations, 'uniform', 0);
            inverseImages = cellfun(@(i) i', images, 'uniform', 0);
            rho = self.repByImages('R', dimension, images, inverseImages);
        end

    end

end
