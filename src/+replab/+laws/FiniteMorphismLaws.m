classdef FiniteMorphismLaws < replab.laws.MorphismLaws

    properties (SetAccess = protected)
        I % (`+replab.FiniteGroup`): Image group
    end

    methods

        function self = FiniteMorphismLaws(morphism)
            self@replab.laws.MorphismLaws(morphism);
            self.I = morphism.image;
        end

% $$$         function law_restrictedSource_S(self, s)
% $$$             sub = self.S.subgroup({s});
% $$$             m = self.morphism.restrictedSource(sub);
% $$$             s1 = sub.sample;
% $$$             self.T.assertEqv(self.morphism.imageElement(s1), m.imageElement(s1));
% $$$         end

% $$$         function isoLaws = laws_toIsomorphism(self)
% $$$             if self.morphism.kernel.isTrivial && ~isa(self.morphism, 'replab.Isomorphism')
% $$$                 isoLaws = self.morphism.toIsomorphism.laws;
% $$$             else
% $$$                 isoLaws = replab.laws.Collection(cell(1, 0));
% $$$             end
% $$$         end
% $$$
% $$$         function law_kernel_(self)
% $$$             gens = self.morphism.kernel.generators;
% $$$             assert(all(cellfun(@(g) self.T.isIdentity(self.morphism.imageElement(g)), gens)));
% $$$         end

        function law_preimageRepresentative_I(self, img)
           s = self.morphism.preimageRepresentative(img);
           img1 = self.morphism.imageElement(s);
           self.T.assertEqv(img, img1);
       end

% $$$        function law_preimagesElement_I(self, img)
% $$$            sset = self.morphism.preimagesElement(img);
% $$$            img1 = self.morphism.imageElement(sset.sample);
% $$$            self.T.assertEqv(img, img1);
% $$$        end
% $$$
% $$$        function law_preimageGroup_(self)
% $$$            K1 = self.morphism.preimageGroup(self.T.trivialSubgroup);
% $$$            K2 = self.morphism.kernel;
% $$$            assertTrue(K1 == K2);
% $$$        end

% $$$        function law_image_(self)
% $$$            gens = cellfun(@(g) self.morphism.imageElement(g), self.morphism.source.generators, 'uniform', 0);
% $$$            img = self.T.subgroup(gens);
% $$$            assertTrue(img == self.morphism.image);
% $$$        end
% $$$
% $$$        function law_imageGroup_S(self, s)
% $$$            GS = self.S.subgroup({s});
% $$$            GT = self.morphism.imageGroup(GS);
% $$$            assert(mod(double(GS.order), double(GT.order)) == 0);
% $$$        end

    end

end
