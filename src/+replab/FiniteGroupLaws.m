classdef FiniteGroupLaws < replab.GroupLaws
    properties (SetAccess = protected)
        I; % index of generator
    end
    methods
        function self = FiniteGroupLaws(T)
            self@replab.GroupLaws(T);
            self.I = replab.domain.intAsDouble(1, T.nGenerators);
        end
    end
    methods
        function law_generatorInverse_I(self, i)
            t = self.T.inverse(self.T.generator(i));
            t1 = self.T.generatorInverse(i);
            self.T.assertEqv(t, t1);
        end
        function law_isTrivial(self)
            self.assert(self.T.isTrivial == (self.T.nGenerators == 0));
        end
        function law_contains_T(self, t)
            self.assert(self.T.contains(t));
        end
        function law_order(self)
            self.assert(self.T.isTrivial == (self.T.order == 1));
        end
        function law_order_elements(self)
            self.assert(self.T.elements.size == self.T.order);
        end
        function elementsLaws = laws_elements(self)
            elementsLaws = replab.IndexedFamilyLaws(self.T.elements);
        end
    end    
end
