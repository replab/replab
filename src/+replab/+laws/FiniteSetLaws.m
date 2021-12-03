classdef FiniteSetLaws < replab.laws.DomainLaws

    properties (SetAccess = protected)
        T % (`+replab.FiniteGroupType`): Finite group type
    end

    methods

        function self = FiniteSetLaws(S)
            self@replab.laws.DomainLaws(S);
            self.T = S.type;
        end

        function yellow(self)
        % Bails out the test for a law if the number of elements > 10000
            if self.S.nElements > 10000
                replab.Laws.skip;
            end
        end

        function red(self)
        % Bails out the test for a law if the number of elements > 100
            if self.S.nElements > 100
                replab.Laws.skip;
            end
        end

        function infrared(self)
        % Bails out the test for a law if the number of elements > 50
            if self.S.nElements > 50
                replab.Laws.skip;
            end
        end

    end

    methods % Laws

        % contains

        function law_set_contains_random_element_S(self, s)
            assertTrue(self.S.contains(s));
        end

        % elements

        function law_elements_(self)
            self.red;
            E = self.S.elements;
            assertTrue(length(E) == self.S.nElements);
            for i = 1:length(E)
                assertTrue(self.S.contains(E{i}));
            end
        end

        % elementsSequence

        function laws = laws_elementsSequence(self)
            self.red;
            laws = self.S.elementsSequence.laws;
        end

        % nElements

        function law_nElements_(self)
            assert(self.S.nElements >= 0);
        end

        % representative

        function law_representative_is_minimal_S(self, s)
            assertTrue(self.S.type.compare(self.S.representative, s) <= 0);
        end

        function law_representative_is_inside_(self)
            assertTrue(self.S.contains(self.S.representative));
        end

        % setProduct

        function laws = laws_setProduct(self)
            laws = self.S.setProduct.laws(self.S);
        end

    end

end
