classdef FiniteGroupTypeLaws < replab.laws.GroupLaws & replab.laws.TotalOrderLaws
% Group axioms

    methods

        function self = FiniteGroupTypeLaws(T)
            self@replab.laws.GroupLaws(T);
            self@replab.laws.TotalOrderLaws(T);
        end

    end

end
