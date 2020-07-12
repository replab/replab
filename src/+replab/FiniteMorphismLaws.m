classdef FiniteMorphismLaws < replab.Laws

    properties (SetAccess = protected)
        morphism % (`+replab.FiniteMorphism`): Morphism tested
        S % `.FiniteGroup`: Source group
        T % `.FiniteGroup`: Target group
    end

    methods

        function self = FiniteMorphismLaws(morphism)
            self.morphism = morphism;
            self.S = morphism.source;
            self.T = morphism.target;
        end

    end

    methods % LAWS

        function law_inverse_S(self, s)
            t = self.morphism.imageElement(s);
            sI = self.S.inverse(s);
            tI1 = self.T.inverse(t);
            tI2 = self.morphism.imageElement(sI);
            self.T.assertEqv(tI1, tI2);
        end

        function law_composition_SS(self, s1, s2)
            s12 = self.S.compose(s1, s2);
            t1 = self.morphism.imageElement(s1);
            t2 = self.morphism.imageElement(s2);
            t12_1 = self.morphism.imageElement(s12);
            t12_2 = self.T.compose(t1, t2);
            self.T.assertEqv(t12_1, t12_2);
        end

        function law_identity_(self)
            self.T.assertEqv(self.T.identity, self.morphism.imageElement(self.S.identity));
        end

    end

end
