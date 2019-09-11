classdef NiceFiniteGroupLaws < replab.FiniteGroupLaws
    methods
        function self = NiceFiniteGroupLaws(T)
            self@replab.FiniteGroupLaws(T);
        end
    end
    methods
        function law_cyclic_subgroup_order_T(self, t)
            if self.T.isIdentity(t)
                sub = self.T.subgroup({});
            else
                sub = self.T.subgroup({t});
            end
            self.assert(sub.order == self.T.elementOrder(t));
        end
    end
end
