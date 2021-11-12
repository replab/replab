classdef FiniteGroupTypeLaws < replab.laws.GroupLaws & replab.laws.TotalOrderLaws
% Group axioms

    methods

        function self = FiniteGroupTypeLaws(T)
            self@replab.laws.GroupLaws(T);
            self@replab.laws.TotalOrderLaws(T);
        end

        function law_is_same_type_as_self_(self)
            self.assert(self.T.isSameTypeAs(self.T));
        end

    end

end
