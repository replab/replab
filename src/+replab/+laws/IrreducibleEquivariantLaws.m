classdef IrreducibleEquivariantLaws < replab.laws.SubEquivariantLaws

    methods

        function self = IrreducibleEquivariantLaws(S)
            self@replab.laws.SubEquivariantLaws(S);
        end

        function law_projectAndFactor_S(self, X)
            [M, err] = self.S.projectAndFactor(X);
            X1 = self.S.reconstruct(M);
            self.assertApproxEqual(X, X1, err);
        end

        function law_projectAndFactorFromParent_S(self, X)
            parentX = self.S.repR.injection * X * self.S.repC.projection;
            [M, err] = self.S.projectAndFactorFromParent(parentX);
            X1 = self.S.reconstruct(M);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
