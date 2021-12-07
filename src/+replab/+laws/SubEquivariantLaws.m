classdef SubEquivariantLaws < replab.laws.EquivariantLaws
% Law checks for a `+replab.SubEquivariant` space

    methods

        function self = SubEquivariantLaws(S)
            self@replab.laws.EquivariantLaws(S);
        end

        function law_projectFromParent_S(self, X)
            parentX = self.S.repR.injection * X * self.S.repC.projection;
            [X1, err] = self.S.projectFromParent(parentX);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
