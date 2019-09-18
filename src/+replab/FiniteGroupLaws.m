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
        function law_order(self)
            self.assert(self.T.isTrivial == (self.T.order == 1));
        end
        function law_order_elements(self)
            self.assert(self.T.elements.size == self.T.order);
        end
        function law_generators(self)
            T = self.T;
            for i = 1:T.nGenerators
                g = T.generator(i);
                ginv = T.generatorInverse(i);
                T.assertEqv(T.identity, T.compose(g, ginv)); % generator consistent with its inverse
                self.assert(T.elements.find(g) > 0); % generator is part of elements
                self.assert(T.elements.find(ginv) > 0); % generator inverse is part of elements
            end
        end
        function elementsLaws = laws_elements(self)
            elementsLaws = replab.IndexedFamilyLaws(self.T.elements);
        end
        function law_decomposition_size(self)
            D = self.T.decomposition.T;
            o = vpi(1);
            for i = 1:length(D)
                o = o * length(D{i});
            end
            self.assert(o == self.T.order);
        end
    end
end
