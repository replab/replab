classdef OrderPreservingFiniteIsomorphismLaws < replab.laws.FiniteIsomorphismLaws

    methods

        function self = OrderPreservingFiniteIsomorphismLaws(morphism)
            self@replab.laws.FiniteIsomorphismLaws(morphism);
        end

        function law_preserves_total_order_SS(self, s1, s2)
            t1 = self.morphism.imageElement(s1);
            t2 = self.morphism.imageElement(s2);
            self.assert(self.S.type.compare(s1, s2) == self.T.type.compare(t1, t2));
        end

    end

end
