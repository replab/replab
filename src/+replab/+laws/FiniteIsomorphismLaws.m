classdef FiniteIsomorphismLaws < replab.laws.IsomorphismLaws & replab.laws.FiniteMorphismLaws

    methods

        function self = FiniteIsomorphismLaws(morphism)
            self@replab.laws.IsomorphismLaws(morphism);
            self@replab.laws.FiniteMorphismLaws(morphism);
        end

    end

    methods % Laws

        function law_preserves_total_order_SS(self, s1, s2)
            if self.morphism.preservesTypeOrder
                t1 = self.morphism.imageElement(s1);
                t2 = self.morphism.imageElement(s2);
                self.assert(self.S.type.compare(s1, s2) == self.T.type.compare(t1, t2));
            end
        end

    end
end
