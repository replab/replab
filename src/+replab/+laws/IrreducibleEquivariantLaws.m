classdef IrreducibleEquivariantLaws < replab.laws.SubEquivariantLaws

    methods

        function self = IrreducibleEquivariantLaws(T)
            self@replab.laws.SubEquivariantLaws(T);
        end

        function law_projectAndFactor_T(self, X)
            [M, err] = self.T.projectAndFactor(X);
            X1 = self.T.reconstruct(M);
            self.assertApproxEqual(X, X1, err);
        end

        function law_projectAndFactorFromParent_T(self, X)
            parentX = self.T.repR.injection * X * self.T.repC.projection;
            [M, err] = self.T.projectAndFactorFromParent(parentX);
            X1 = self.T.reconstruct(M);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
