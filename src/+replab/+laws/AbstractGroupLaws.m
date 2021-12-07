classdef AbstractGroupLaws < replab.laws.FiniteGroupLaws
% Law checks for abstract groups: the operations below are pretty expensive

    methods

        function self = AbstractGroupLaws(S)
            self@replab.laws.FiniteGroupLaws(S);
        end

    end

    methods

        function law_simplify_SS(self, s1, s2)
            s = self.S.compose(s1, s2);
            s1 = self.S.simplify(s);
            self.S.assertEqv(s, s1);
        end

        function law_withGeneratorNames_S(self, s)
        % nothing, disable this law
        end

    end

end
