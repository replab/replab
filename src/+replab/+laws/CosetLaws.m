classdef CosetLaws < replab.laws.FiniteSetLaws

    methods

        function self = CosetLaws(T)
            self@replab.laws.FiniteSetLaws(T);
        end

    end

    methods

        function law_factorize_coset_representative_(self)
            [l, r] = self.T.factorizeShortRepresentativeLetters;
            r1 = self.T.parent.imageLetters(l);
            self.T.parent.assertEqv(r, r1);
            assertTrue(self.T.contains(r));
        end

    end

end
