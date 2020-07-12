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
% The `.niceMonomorphismImage` method of this group should be valid for all elements of the parent as well;
% then the `.contains` method will be valid for all elements of the parent group.

    properties (SetAccess = protected)

        niceGroup % (`+replab.PermutationGroup`): Permutation group isomorphic to this finite group
        niceMorphism % (`+replab.FiniteIsomorphism`): Isomorphism from this group to `.niceGroup`

    end


    %% Abstract

    methods

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
            chain = replab.bsgs.Chain.make(n, generators);
            P = replab.PermutationGroup(n, generators, chain.order, self.niceGroup.parent, chain);
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

        function order = computeOrder(self)
            order = self.niceGroup.chain.order;
        end

        function m = computeNiceInverseMonomorphism(self)
            m = replab.mrp.PermMorphism.byImages(self.niceGroup, self, self.generators);
        end

        function g = computeNiceGroup(self)
            imgId = self.niceMonomorphismImage(self.identity);
            n = length(imgId);
            g = replab.PermutationGroup(n, self.niceGenerators, self.cachedOrEmpty('order'));
        end

        function dec = computeDecomposition(self)
            dec = replab.FiniteGroupDecomposition(self, self.niceInverseMonomorphism.chain.imagesDecomposition);
        end

    end

    methods

        function c = niceGroup(self)
        % Returns the image of this group as a permutation group through the nice monomorphismx
        %
        % Returns:
        %   `+replab.PermutationGroup`: Permutation group
            c = self.cached('niceGroup', @() self.computeNiceGroup);
        end

        function m = niceInverseMonomorphism(self)
        % Returns the monomorphism from the permutation representation to the original group
            m = self.cached('niceInverseMonomorphism', @() self.computeNiceInverseMonomorphism);
        end

    end

    methods % Computed properties


    end

    methods

        function f = leftConjugateMorphism(self, by)
        % Returns the morphism that corresponds to left conjugation by an element
        %
        % Args:
        %   by (element of `parent`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.Morphism`: Conjugation morphism
            generatorImages = cellfun(@(g) self.parent.leftConjugate(by, g), self.generators, 'uniform', 0);
            target = self.parent.subgroup(newGenerators, self.order);
            f = self.morphismByImages(self, target, generatorImages);
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
        %   >>> S4 = replab.SymmetricGroup(4);
        %   >>> S4.order
        %       24
        %   >>> D = S4.derivedSubgroup;
        %   >>> D == replab.AlternatingGroup(4)
        %       1
        % Returns:
        %   `+replab.NiceFiniteGroup`: The derived subgroup
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            sub = self.niceMonomorphismGroupPreimage(self.niceGroup.derivedSubgroup);
        end

        function c = rightCosetsOf(self, subgroup)
            assert(~isa(self, 'replab.PermutationGroup')); % handled in subclass
            c = replab.nfg.RightCosets(self, subgroup);
        end

        function c = mldivide(self, supergroup)
            c = supergroup.rightCosetsOf(self);
        end

        function c = leftCosetsOf(self, subgroup)
            assert(~isa(self, 'replab.PermutationGroup')); % handled in subclass
            c = replab.nfg.LeftCosets(self, subgroup);
        end

        function c = mrdivide(self, subgroup)
            c = self.leftCosetsOf(subgroup);
        end

        function sub = centralizer(self, rhs)
        % Returns the centralizer of the given element in this group
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
        %   rhs (`.NiceFiniteGroup` or group element): Element to compute the centralizer of
        %
        % Returns:
        %   `.NiceFiniteGroup`: The centralizer
        %
        % See `.centralizerElement` and `.centralizerGroup` which it calls.
            if isa(rhs, 'replab.NiceFiniteGroup')
                sub = self.centralizerGroup(rhs);
            else
                sub = self.centralizerElement(rhs);
            end
        end

        function sub = centralizerElement(self, g)
        % Returns the centralizer of a group in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerElement([2 3 1 4]);
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Args:
        %   g (group element): Element to compute the centralizer of
        %
        % Returns:
        %   `.NiceFiniteGroup`: The centralizer
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            sub = self.niceMonomorphismGroupPreimage(self.niceGroup.centralizerElement(self.niceMonomorphismImage(g)));
        end

        function sub = centralizerGroup(self, G)
        % Returns the centralizer of a group in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerGroup(G.subgroup({[2 3 1 4]}));
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerGroup(G.subgroup({[2 3 1 4] [2 1 3 4]}));
        %   >>> C == replab.S(4).trivialGroup
        %     1
        %
        % Args:
        %   G (`.NiceFiniteGroup`): Permutation group with the same parent as this group
        %
        % Returns:
        %   `.NiceFiniteGroup`: The centralizer
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            sub = self.niceMonomorphismGroupPreimage(self.niceGroup.centralizerGroup(self.niceMonomorphismGroupImage(G)));
        end

        function res = closure(self, obj)
        % Computes the group generated by the elements of this group and another given element
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
        %   rhs (`.NiceFiniteGroup` or group element): Element(s) to add
        %
        % Returns:
        %   `.NiceFiniteGroup`: The closure
            if isa(obj, 'replab.NiceFiniteGroup')
                res = self.closureGroup(obj);
            else
                res = self.closureElement(obj);
            end
        end

        function res = closureElement(self, g)
        % Computes the group generated by the elements of this group and another given element
        %
        % Args:
        %   g (element of `.parent`): Element to add
        %
        % Returns:
        %   `.NiceFiniteGroup`: The closure
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            res = self.parent.niceMonomorphismGroupPreimage(self.niceGroup.closureElement(self.niceMonomorphismImage(g)));
        end

        function res = closureGroup(self, G)
        % Computes the group generated by the elements of this group and another given element
        %
        % Args:
        %   G (`.NiceFiniteGroup` with same parent): Elements to add
        %
        % Returns:
        %   `.NiceFiniteGroup`: The closure
            assert(~isa(self, 'replab.PermutationGroup')); % is handled in subclass
            res = self.parent.niceMonomorphismGroupPreimage(self.niceGroup.closureGroup(self.niceMonomorphismImage(g)));
        end


        function res = isCyclic(self)
        % Returns whether this group is a cyclic group
            assert(~isa(self, 'replab.PermutationGroup'));
            if isprime(self.order)
                res = true;
            else
                res = self.niceGroup.isCyclic;
            end
        end

        function grp = trivialSubgroup(self)
        % Returns the trivial subgroup of this group
        %
        % Returns:
        %   `+replab.NiceFiniteGroup`: The trivial subgroup
            grp = self.subgroup({}, vpi(1));
        end

        function c = conjugacyClasses(self)
        % Returns the conjugacy classes of this group
        %
        % Returns:
        %   cell(1, \*) of `+replab.ConjugacyClass`: Array of conjugacy classes
            c = replab.ConjugacyClass.computeAll(self);
        end

        function o = elementOrder(self, g)
        % Returns the order of a group element
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   vpi: The order of ``g``, i.e. the smallest ``o`` such that ``g^o == identity``
            o = self.niceGroup.elementOrder(self.niceMonomorphismImage(g));
        end

        %% Methods enabled by the BSGS algorithms

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

        function g = sample(self)
            [~, g] = self.niceInverseMonomorphism.chain.sample;
        end

    end

end
