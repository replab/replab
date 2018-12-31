classdef FiniteGroupLaws < replab.FinitelyGeneratedGroupLaws
    methods
        function self = FiniteGroupLaws(T)
            self = self@replab.FinitelyGeneratedGroupLaws(T);
        end
    end
    methods
        function law_contains_T(self, t)
            self.assert(self.T.contains(t));
        end
        function law_order(self)
            self.assert(self.T.isTrivial == (self.T.order == 1));
        end
        function law_order_elements(self)
            self.assert(self.T.elements.size == self.T.order);
        end
        function enumeratorLaws = laws_enumerator(self)
            enumeratorLaws = replab.EnumeratorLaws(self.T.elements);
        end
    end    
end
