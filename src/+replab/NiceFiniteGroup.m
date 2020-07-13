classdef NiceFiniteGroup < replab.FiniteGroup
% A nice finite group is a finite group equipped with an injective homomorphism into a permutation group
%
% Any class that subclasses `.NiceFiniteGroup` must implements a method `.niceImage` that returns a permutation
% corresponding to a group element.
%
% Then, an isomorphism is established between the present `.NiceFiniteGroup` and a permutation group; as
% permutation groups can be handled by efficient BSGS algorithms, the requested computations can be
% translated back and forth between this group and a permutation group.
%
% In particular, the decomposition of the finite group in a product of sets (`.FiniteGroupDecomposition`),
% the enumeration of elements using a `.IndexedFamily`, the construction of subgroups is all handled
% by permutation group algorithms.
%
% Note that the `.niceImage` method of this group should be valid for all elements of the parent as well,
% and that groups having the same type (as verified by `.hasSameTypeAs`) should return the same images
% under `.niceImage`.

    methods % Nice monomorphism support

        function p = niceImage(self, g)
        % Returns a permutation image of the given group element
        %
        % Args:
        %   g (element): Group element to represent as a permutation
        %
        % Returns:
        %   permutation: Permutation image of ``g``
            error('Abstract');
        end

        function G = niceGroup(self)
        % Returns the permutation group isomorphic to this finite group
        %
        % Returns:
        %   `+replab.PermutationGroup`: The image of the isomorphism
            G = self.cached('niceGroup', @() self.niceMorphism.image);
        end

        function f = niceMorphism(self)
        % Returns the isomorphism from this group to `.niceGroup`
            f = self.cached('niceMorphism', @() replab.nfg.NiceFiniteGroupIsomorphism(self));
        end

    end

    methods

        function order = computeOrder(self)
            order = self.niceGroup.order;
        end

        function dec = computeDecomposition(self)
        % TODO
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
