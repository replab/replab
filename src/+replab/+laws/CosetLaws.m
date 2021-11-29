classdef CosetLaws < replab.laws.FiniteSetLaws

    methods

        function self = CosetLaws(S)
            self@replab.laws.FiniteSetLaws(S);
        end

    end

    methods

        function law_factorize_coset_representative_(self)
            [l, r] = self.S.factorizeShortRepresentativeLetters;
            r1 = self.S.group.imageLetters(l);
            self.S.type.assertEqv(r, r1);
            assertTrue(self.S.contains(r));
        end

    end

end
