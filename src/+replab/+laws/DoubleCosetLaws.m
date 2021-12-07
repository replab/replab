classdef DoubleCosetLaws < replab.laws.FiniteSetLaws

    properties (SetAccess = protected)
        L % (`+replab.FiniteGroup`): Left subgroup
        R % (`+replab.FiniteGroup`): Right subgroup
    end

    methods

        function self = DoubleCosetLaws(S)
            self@replab.laws.FiniteSetLaws(S);
            self.L = S.leftSubgroup;
            self.R = S.rightSubgroup;
        end

    end

    methods

        function law_same_double_coset_LR(self, l, r)
            assert(self.S.contains(self.S.group.composeAll({l, self.S.sample, r})));
        end

    end

end
