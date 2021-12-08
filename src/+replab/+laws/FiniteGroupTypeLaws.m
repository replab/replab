classdef FiniteGroupTypeLaws < replab.laws.GroupLaws & replab.laws.TotalOrderLaws
% Group axioms

    methods

        function self = FiniteGroupTypeLaws(S)
            self@replab.laws.GroupLaws(S);
            self@replab.laws.TotalOrderLaws(S);
        end

        function law_is_same_type_as_self_(self)
            self.assert(self.S.isSameTypeAs(self.S));
        end

    end

end
