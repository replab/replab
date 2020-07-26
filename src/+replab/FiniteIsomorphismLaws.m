classdef FiniteIsomorphismLaws < replab.FiniteMorphismLaws

    methods

        function self = FiniteIsomorphismLaws(morphism)
            self@replab.FiniteMorphismLaws(morphism);
        end

        function law_roundtrip_source_S(self, s)
            t = self.morphism.imageElement(s);
            s1 = self.morphism.preimageElement(t);
            self.S.assertEqv(s, s1);
        end

        function law_roundtrip_target_T(self, t)
            s = self.morphism.preimageElement(t);
            t1 = self.morphism.imageElement(s);
            self.T.assertEqv(t, t1);
        end

    end

end
