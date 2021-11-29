classdef FiniteGroupLaws < replab.laws.GroupLaws & replab.laws.FiniteSetLaws
% Law checks for finite groups: the operations below should be cheap even for large groups

    methods

        function self = FiniteGroupLaws(S)
            self@replab.laws.GroupLaws(S);
            self@replab.laws.FiniteSetLaws(S);
        end

    end

    methods

% $$$         function law_generatorInverse_(self)
% $$$         % Checks generator inverses
% $$$             for i = 1:self.S.nGenerators
% $$$                 t = self.S.inverse(self.S.generator(i));
% $$$                 t1 = self.S.generatorInverse(i);
% $$$                 self.S.assertEqv(t, t1);
% $$$             end
% $$$         end
% $$$
% $$$         function law_isTrivial_(self)
% $$$         % Checks that a group is trivial iff it has no generators
% $$$             self.assert(self.S.isTrivial == (self.S.nGenerators == 0));
% $$$         end
% $$$
% $$$         function law_withGeneratorNames_(self)
% $$$         % Checks that renaming generators works properly
% $$$             G = self.S.withGeneratorNames(replab.fp.defaultGeneratorNames(self.S.nGenerators));
% $$$             assert(G.nGenerators == self.S.nGenerators);
% $$$         end
% $$$
% $$$         function law_order_(self)
% $$$         % Checks that a group is trivial iff its order is 1
% $$$             self.assert(self.S.isTrivial == (self.S.order == 1));
% $$$         end
% $$$
% $$$         function law_order_elements_(self)
% $$$         % Checks that the number of elements corresponds to the group order
% $$$             self.assert(self.S.elementsSequence.nElements == self.S.order);
% $$$         end
% $$$
% $$$         function law_generators_(self)
% $$$         % Performs various safety checks on the group generators
% $$$             T = self.S;
% $$$             for i = 1:T.nGenerators
% $$$                 g = T.generator(i);
% $$$                 ginv = T.generatorInverse(i);
% $$$                 T.assertEqv(T.identity, T.compose(g, ginv)); % generator consistent with its inverse
% $$$                 self.assert(T.elementsSequence.find(g) > 0); % generator is part of elements
% $$$                 self.assert(T.elementsSequence.find(ginv) > 0); % generator inverse is part of elements
% $$$             end
% $$$         end
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
% $$$         function law_cyclic_subgroup_order_S(self, t)
% $$$             self.thisIsSlow;
% $$$             if self.S.isIdentity(t)
% $$$                 sub = self.S.subgroup({});
% $$$             else
% $$$                 sub = self.S.subgroup({t});
% $$$             end
% $$$             self.assert(sub.order == self.S.elementOrder(t));
% $$$         end
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

end
