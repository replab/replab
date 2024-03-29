classdef FiniteGroupLaws < replab.laws.GroupLaws & replab.laws.FiniteSetLaws
% Law checks for finite groups: the operations below should be cheap even for large groups

    methods

        function self = FiniteGroupLaws(S)
            self@replab.laws.GroupLaws(S);
            self@replab.laws.FiniteSetLaws(S);
        end

    end

    methods % Laws: Implementations

        function law_eq_to_itself_(self)
            self.assert(self.S == self.S);
        end

        function law_isequal_to_itself_(self)
            self.assert(isequal(self.S, self.S));
        end

    end

    methods % Laws: Conjugacy classes and character table

        % skipping characterTable, complexCharacterTable, realCharacterTable

        function law_conjugacyClass_S(self, s)
            cl = self.S.conjugacyClass(s);
            self.assert(cl.contains(s));
        end

        function laws = laws_conjugacyClasses(self)
            self.yellow;
            laws = self.S.conjugacyClasses.laws;
        end

        function law_findLeftConjugations_SS(self, b, s)
            t = self.S.leftConjugate(b, s);
            b1 = self.S.findLeftConjugations(s, t);
            assert(b1.contains(b));
        end

    end

    methods % Laws: Group properites

        function law_abelianInvariants_(self)
            self.assert(prod(self.S.abelianInvariants) == self.S.order / self.S.derivedSubgroup.order);
        end

        function law_exponent_S(self, s)
            self.yellow;
            e = self.S.exponent;
            if e > 2^52
                return
            end
            e = double(e);
            self.assert(self.S.isIdentity(self.S.composeN(s, e)));
        end

        function law_isCommutative_SS(self, s, t)
            if self.S.isCommutative
                self.S.assertEqv(self.S.compose(s, t), self.S.compose(t, s));
            end
        end

        function law_isCyclic_(self)
            self.yellow;
            if self.S.isCyclic
                self.assert(self.S.order == self.S.exponent);
            end
        end

        function law_isSimple_(self)
            self.yellow;
            if self.S.isSimple
                self.assert(self.S.derivedSubgroup.isTrivial);
            end
        end

        function law_isTrivial_(self)
            self.assert(self.S.isTrivial == (self.S.order == 1));
        end

        function law_knownOrder_(self)
            self.S.order;
            self.assert(self.S.knownOrder);
        end

        function law_order_(self)
            self.yellow;
            o = self.S.order;
            e = self.S.exponent;
            % exponent divides the order
            self.assert(mod(o, e) == 0);
        end

        function law_presentation_(self)
            self.red;
            G = replab.AbstractGroup.fromPresentation(self.S.presentation);
            self.assert(G.order == self.S.order);
        end

    end

    methods % Laws: Generators-related methods

        function law_factorizeFlat_left_coset_S(self, s)
            subgroup = self.S.randomSubgroup;
            coset = subgroup.leftCoset(s, 'group', self.S);
            l = self.S.factorizeFlat(coset);
            r = self.S.imageFlat(l);
            self.assert(coset.contains(r));
        end

        function law_factorizeFlat_right_coset_S(self, s)
            subgroup = self.S.randomSubgroup;
            coset = subgroup.rightCoset(s, 'group', self.S);
            l = self.S.factorizeFlat(coset);
            r = self.S.imageFlat(l);
            self.assert(coset.contains(r));
        end

        function law_factorizeFlat_imageFlat_element_S(self, s)
            self.S.assertEqv(self.S.imageFlat(self.S.factorizeFlat(s)), s);
        end

        function law_factorizeWord_imageWord_element_S(self, s)
            self.S.assertEqv(self.S.imageWord(self.S.factorizeWord(s)), s);
        end

        function law_knownRelators_(self)
            self.red;
            self.S.relatorsFlat;
            self.assert(self.S.knownRelators);
        end

        function law_flatToWord_wordToFlat_S(self, s)
            w = self.S.factorizeWord(s);
            f = self.S.wordToFlat(w);
            w1 = self.S.flatToWord(f);
            self.S.assertEqv(self.S.imageWord(w1), s);
        end

        function law_generatorInverse_(self)
        % Checks generator inverses
            for i = 1:self.S.nGenerators
                t = self.S.inverse(self.S.generator(i));
                t1 = self.S.generatorInverse(i);
                self.S.assertEqv(t, t1);
            end
        end

        function law_generatorNames_(self)
            names = self.S.generatorNames;
            assert(length(names) == self.S.nGenerators);
            for i = 1:length(names)
                self.S.assertEqv(self.S.imageWord(names{i}), self.S.generator(i));
            end
        end

        function law_generators_(self)
        % Performs various safety checks on the group generators
            T = self.S;
            for i = 1:T.nGenerators
                g = T.generator(i);
                ginv = T.generatorInverse(i);
                T.assertEqv(T.identity, T.compose(g, ginv)); % generator consistent with its inverse
                self.assert(T.elementsSequence.find(g) > 0); % generator is part of elements
                self.assert(T.elementsSequence.find(ginv) > 0); % generator inverse is part of elements
            end
        end

        function law_relatorsFlat_(self)
            self.red;
            R = self.S.relatorsFlat;
            for i = 1:length(R)
                self.assert(self.S.isIdentity(self.S.imageFlat(R{i})));
            end
        end

        function law_relatorsWord_(self)
            self.red;
            R = self.S.relatorsWord;
            for i = 1:length(R)
                self.assert(self.S.isIdentity(self.S.imageWord(R{i})));
            end
        end

    end

    methods % Laws: Elements

        function law_elementOrder_S(self, s)
            self.yellow;
            e = self.S.elementOrder(s);
            self.assert(self.S.isIdentity(self.S.composeN(s, e)));
        end

        function law_smallGeneratingSet_(self)
            self.yellow;
            G = self.S.subgroupWithGenerators(self.S.smallGeneratingSet);
            self.assert(self.S.order == G.order);
            assert(G.nGenerators <= self.S.nGenerators);
        end

    end


    methods % Laws: Construction of groups

        function law_closure_S(self, s)
            G1 = self.S.closure(s);
            G2 = self.S.closure(self.S);
            G3 = self.S.closure(s, self.S);
            G4 = self.S.closure(s, self.S, s);
            self.assert(G1 == self.S);
            self.assert(G2 == self.S);
            self.assert(G3 == self.S);
            self.assert(G4 == self.S);
        end

        function law_leftConjugateGroup_S(self, s)
            self.assert(self.S == self.S.leftConjugateGroup(s));
        end

        function law_normalClosure_(self)
            self.assert(self.S.derivedSubgroup == self.S.normalClosure(self.S.derivedSubgroup));
        end

        function law_withGeneratorNames_S(self, s)
            if self.S.isIdentity(s)
                return
            end
            G = self.S.subgroup({s});
            G1 = G.withGeneratorNames({'x'});
            self.S.assertEqv(G1.imageWord('x'), s);
        end

    end

    methods % Laws: Subgroups

        function law_center_(self)
            self.assert(self.S.center.isCommutative);
        end

        function law_centralizer_syntax_S(self, s)
            G1 = self.S.centralizer(s);
            G2 = self.S.centralizer(self.S.subgroup({s}));
            self.assert(G1 == G2);
        end

        function law_centralizer_identity_(self)
            self.assert(self.S.centralizer(self.S.identity) == self.S);
        end

        % derivedSeries has a simple implementation, no test

        function law_derivedSubgroup_(self)
            self.assert(self.S.derivedSubgroup.isNormalSubgroupOf(self.S));
        end

        function laws = laws_subgroups_lattice(self)
            self.yellow;
            % checks closure, intersection, randomProperSubgroup
            laws = replab.domain.Subgroups(self.S).laws;
        end

        % subgroup

        function law_cyclic_subgroup_order_S(self, t)
            self.yellow;
            if self.S.isIdentity(t)
                sub = self.S.subgroup({});
            else
                sub = self.S.subgroup({t});
            end
            self.assert(sub.order == self.S.elementOrder(t));
        end

        function law_subgroupWithGenerators_S(self, s)
            if self.S.isIdentity(s)
                return
            end
            G1 = self.S.subgroup({s});
            G2 = self.S.subgroupWithGenerators({s});
            assert(G1 == G2);
        end

        function law_trivialSubgroup_(self)
            self.assert(self.S.trivialSubgroup.isTrivial);
        end

    end

    methods % Laws: Cosets

        function law_doubleCoset_SS(self, g, h)
            self.red;
            S = self.S;
            T = S.trivialSubgroup;
            C1 = S.doubleCoset(g, T);
            C2 = T.doubleCoset(g, S);
            C3 = S.doubleCoset(g, S);
            self.assert(C1.contains(h));
            self.assert(C2.contains(h));
            self.assert(C3.contains(h));
        end

        function laws = laws_doubleCoset(self)
            self.red;
            L = self.S.randomSubgroup;
            R = self.S.randomSubgroup;
            el = self.S.sample;
            dc = L.doubleCoset(el, R, 'group', self.S);
            laws = dc.laws;
        end

        function law_doubleCosets_(self)
            self.red;
            C = self.S.doubleCosets(self.S, self.S);
            self.assert(C.nElements == 1);
        end

        function laws = laws_doubleCosets(self)
            self.red;
            laws = self.S.doubleCosets(self.S.randomSubgroup, self.S.randomSubgroup).laws;
        end

        function law_isNormalizedBy_S(self, g)
            self.assert(self.S.isNormalizedBy(g));
        end

        function laws = laws_leftCoset(self)
            self.red;
            el = self.S.sample;
            G = self.S.randomSubgroup;
            lc = G.leftCoset(el, 'group', self.S);
            laws = lc.laws;
        end

        function laws = laws_leftCosets(self)
            self.red;
            laws = self.S.leftCosets(self.S.randomSubgroup).laws;
        end

        % we skip mldivide and mrdivide as they relate to syntax

        function laws = laws_normalCoset(self)
            self.red;
            G = self.S.randomSubgroup;
            nc = G.normalCoset(G.identity, 'group', self.S);
            assert(isa(nc, 'replab.NormalCoset'));
            laws = nc.laws;
        end

        function laws = laws_normalCosets(self)
            self.red;
            laws = self.S.normalCosets(self.S.derivedSubgroup).laws;
        end

        function laws = laws_rightCoset(self)
            self.red;
            el = self.S.sample;
            G = self.S.randomSubgroup;
            rc = G.rightCoset(el, 'group', self.S);
            laws = rc.laws;
        end

        function laws = laws_rightCosets(self)
            self.red;
            laws = self.S.rightCosets(self.S.randomSubgroup).laws;
        end

    end

    methods % Relations to other groups

        function law_derived_subgroup_is_normal_subgroup_(self)
            self.assert(self.S.derivedSubgroup.isNormalSubgroupOf(self.S));
        end

        function law_isSubgroupOf_(self)
            self.assert(self.S.randomSubgroup.isSubgroupOf(self.S));
        end

    end

    methods % Isomorphic groups

        function law_abstractGroup_(self)
            self.red;
            T = self.S.abstractGroup;
            iso = self.S.abstractIsomorphism;
            self.assert(self.S.nGenerators == T.nGenerators);
            for i = 1:self.S.nGenerators
                T.assertEqv(T.generator(i), iso.imageElement(self.S.generator(i)));
            end
        end

        function law_permutationGroup_(self)
            self.yellow;
            P = self.S.permutationGroup;
            iso = self.S.permutationIsomorphism;
            self.assert(self.S.nGenerators == P.nGenerators);
            for i = 1:self.S.nGenerators
                P.assertEqv(P.generator(i), iso.imageElement(self.S.generator(i)));
            end
        end

    end

    methods % Morphisms

        function laws = laws_abstractIsomorphism(self)
            self.red;
            laws = self.S.abstractIsomorphism.laws;
        end

        function laws = laws_conjugatingAutomorphism(self)
            self.red;
            laws = self.S.conjugatingAutomorphism(self.S.sample).laws;
        end

        function law_findIsomorphism_to_trivial_group_(self)
            self.infrared;
            tg = replab.PermutationGroup.trivial(2);
            self.assert(isempty(self.S.findIsomorphism(tg)) ~= self.S.isTrivial);
        end

        function law_findIsomorphisms_(self)
            self.infrared;
            nInner = double(self.S.order / self.S.center.order);
            self.assert(length(self.S.findIsomorphisms(self.S)) >= nInner);
        end

        % findMorphisms

        function law_findMorphisms_(self)
            self.infrared;
            tg = replab.PermutationGroup.trivial(2);
            % there is only one morphism to the trivial group
            self.assert(length(self.S.findMorphisms(tg)) == 1);
        end

        function law_isMorphismByImages_(self)
            self.red;
            self.assert(self.S.isMorphismByImages(self.S, 'preimages', self.S.generators, 'images', self.S.generators));
        end

        % TODO isomorphismByFunction, isomorphismByImages, morphismByImages

        function law_orderPreservingPermutationIsomorphism_(self)
            self.red;
            iso = self.S.orderPreservingPermutationIsomorphism;
            self.assert(iso.preservesTypeOrder);
        end

        function laws = laws_orderPreservingPermutationIsomorphism(self)
            self.red;
            laws = self.S.orderPreservingPermutationIsomorphism.laws;
        end

        function laws = laws_permutationIsomorphism(self)
            self.red;
            laws = self.S.permutationIsomorphism.laws;
        end

    end

% $$$
% $$$
% $$$         function law_withGeneratorNames_(self)
% $$$         % Checks that renaming generators works properly
% $$$             G = self.S.withGeneratorNames(replab.fp.defaultGeneratorNames(self.S.nGenerators));
% $$$             assert(G.nGenerators == self.S.nGenerators);
% $$$         end
% $$$
% $$$         function law_order_elements_(self)
% $$$         % Checks that the number of elements corresponds to the group order
% $$$             self.assert(self.S.elementsSequence.nElements == self.S.order);
% $$$         end
% $$$
% $$$
% $$$         function elementsLaws = laws_elementsSequence(self)
% $$$         % Tests the group elements as an indexed family
% $$$             elementsLaws = self.S.elementsSequence.laws;
% $$$         end
% $$$
% $$$         function law_setProduct_size_(self)
% $$$         % Checks that the cartesian product set decomposition has the correct size
% $$$             D = self.S.setProduct.sets;
% $$$             o = vpi(1);
% $$$             for i = 1:length(D)
% $$$                 o = o * length(D{i});
% $$$             end
% $$$             self.assert(o == self.S.order);
% $$$         end
% $$$
% $$$         function law_contains_S(self, t)
% $$$         % Checks element membership (trivial case)
% $$$             self.assert(self.S.contains(t));
% $$$         end
% $$$
% $$$
% $$$         function law_factorizeFlat_imageFlat_S(self, t)
% $$$             l = self.S.factorizeFlat(t);
% $$$             t1 = self.S.imageFlat(l);
% $$$             self.S.assertEqv(t, t1);
% $$$         end
% $$$
% $$$         function law_factorizeWord_imageWord_S(self, t)
% $$$             w = self.S.factorizeWord(t);
% $$$             t1 = self.S.imageWord(w);
% $$$             self.S.assertEqv(t, t1);
% $$$         end
% $$$
% $$$         function law_relators_are_satisfied_(self)
% $$$             if self.S.knownRelators
% $$$                 rels = self.S.relatorsFlat;
% $$$                 for i = 1:length(rels)
% $$$                     g = self.S.imageFlat(rels{i});
% $$$                     assert(self.S.isIdentity(g));
% $$$                 end
% $$$             end
% $$$         end

end
