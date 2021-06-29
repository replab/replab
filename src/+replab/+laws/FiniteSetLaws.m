classdef FiniteSetLaws < replab.laws.DomainLaws

    methods

        function self = FiniteSetLaws(T)
            self@replab.laws.DomainLaws(T);
        end

    end

    methods % Laws

        function law_set_contains_representative_(self)
            assertTrue(self.T.contains(self.T.representative));
        end

        function law_set_contains_random_element_T(self, t)
            assertTrue(self.T.contains(t));
        end

    end

end
