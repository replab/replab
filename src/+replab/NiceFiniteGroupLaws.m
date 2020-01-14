classdef NiceFiniteGroupLaws < replab.FiniteGroupLaws
% Law checks for nice finite groups
    methods

        function self = NiceFiniteGroupLaws(T)
            self@replab.FiniteGroupLaws(T);
        end

    end

    methods

        function law_contains_T(self, t)
        % Checks element membership (trivial case)
            self.assert(self.T.contains(t));
        end
% $$$         function law_cyclic_subgroup_order_T(self, t)
% $$$             if self.T.isIdentity(t)
% $$$                 sub = self.T.subgroup({});
% $$$             else
% $$$                 sub = self.T.subgroup({t});
% $$$             end
% $$$             self.assert(sub.order == self.T.elementOrder(t));
% $$$         end
    end
end
