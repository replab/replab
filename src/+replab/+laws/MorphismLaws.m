classdef MorphismLaws < replab.Laws

    properties (SetAccess = protected)
        morphism % (`+replab.Morphism`): Morphism tested
        S % (`+replab.Group`): Source group
        T % (`+replab.Group`): Target group
    end

    methods

        function self = MorphismLaws(morphism)
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

        function L = laws_torusMap(self)
            if ~isempty(self.morphism.torusMap)
                L = replab.laws.TorusMorphismLaws(self.morphism);
            else
                L = replab.Laws.empty;
            end
        end

    end

end
