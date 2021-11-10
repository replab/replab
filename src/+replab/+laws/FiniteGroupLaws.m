classdef FiniteGroupLaws < replab.laws.GroupLaws
% Law checks for finite groups: the operations below are pretty expensive

    methods

        function self = FiniteGroupLaws(T)
            self@replab.laws.GroupLaws(T);
        end

    end

    methods

        function law_generatorInverse_(self)
        % Checks generator inverses
            for i = 1:self.T.nGenerators
                t = self.T.inverse(self.T.generator(i));
                t1 = self.T.generatorInverse(i);
                self.T.assertEqv(t, t1);
            end
        end

        function law_isTrivial_(self)
        % Checks that a group is trivial iff it has no generators
            self.assert(self.T.isTrivial == (self.T.nGenerators == 0));
        end

        function law_withGeneratorNames_(self)
        % Checks that renaming generators works properly
            G = self.T.withGeneratorNames(replab.fp.defaultGeneratorNames(self.T.nGenerators));
            assert(G.nGenerators == self.T.nGenerators);
        end

        function law_order_(self)
        % Checks that a group is trivial iff its order is 1
            self.assert(self.T.isTrivial == (self.T.order == 1));
        end

        function law_order_elements_(self)
        % Checks that the number of elements corresponds to the group order
            self.assert(self.T.elementsSequence.nElements == self.T.order);
        end

        function law_generators_(self)
        % Performs various safety checks on the group generators
            T = self.T;
            for i = 1:T.nGenerators
                g = T.generator(i);
                ginv = T.generatorInverse(i);
                T.assertEqv(T.identity, T.compose(g, ginv)); % generator consistent with its inverse
                self.assert(T.elementsSequence.find(g) > 0); % generator is part of elements
                self.assert(T.elementsSequence.find(ginv) > 0); % generator inverse is part of elements
            end
        end

        function elementsLaws = laws_elementsSequence(self)
        % Tests the group elements as an indexed family
            elementsLaws = self.T.elementsSequence.laws;
        end

        function law_setProduct_size_(self)
        % Checks that the cartesian product set decomposition has the correct size
            D = self.T.setProduct.sets;
            o = vpi(1);
            for i = 1:length(D)
                o = o * length(D{i});
            end
            self.assert(o == self.T.order);
        end

        function law_contains_T(self, t)
        % Checks element membership (trivial case)
            self.assert(self.T.contains(t));
        end

        function law_cyclic_subgroup_order_T(self, t)
            self.thisIsSlow;
            if self.T.isIdentity(t)
                sub = self.T.subgroup({});
            else
                sub = self.T.subgroup({t});
            end
            self.assert(sub.order == self.T.elementOrder(t));
        end

        function law_factorizeFlat_imageFlat_T(self, t)
            l = self.T.factorizeFlat(t);
            t1 = self.T.imageFlat(l);
            self.T.assertEqv(t, t1);
        end

        function law_factorizeWord_imageWord_T(self, t)
            w = self.T.factorizeWord(t);
            t1 = self.T.imageWord(w);
            self.T.assertEqv(t, t1);
        end

        function law_relators_are_satisfied_(self)
            if self.T.knownRelators
                rels = self.T.relatorsFlat;
                for i = 1:length(rels)
                    g = self.T.imageFlat(rels{i});
                    assert(self.T.isIdentity(g));
                end
            end
        end

    end

end
