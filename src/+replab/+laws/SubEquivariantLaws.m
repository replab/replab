classdef SubEquivariantLaws < replab.laws.EquivariantLaws
% Law checks for a `+replab.SubEquivariant` space

    methods

        function self = SubEquivariantLaws(T)
            self@replab.laws.EquivariantLaws(T);
        end

        function law_projectFromParent_T(self, X)
            parentX = self.T.repR.injection * X * self.T.repC.projection;
            [X1, err] = self.T.projectFromParent(parentX);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
