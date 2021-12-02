classdef CosetLaws < replab.laws.FiniteSetLaws

    properties (SetAccess = protected)
        G % (`+replab.FiniteGroup`): Subgroup
    end

    methods

        function self = CosetLaws(S)
            self@replab.laws.FiniteSetLaws(S);
            self.G = S.subgroup;
        end

    end

    methods

        function law_factorize_coset_representative_(self)
            l = self.S.group.factorizeFlat(self.S);
            r = self.S.group.imageFlat(l);
            assertTrue(self.S.contains(r));
        end

    end

end
