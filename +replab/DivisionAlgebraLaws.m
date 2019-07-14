classdef DivisionAlgebraLaws < replab.MonoidLaws
    methods
        function self = DivisionAlgebraLaws(T)
            self@replab.MonoidLaws(T);
        end
        function law_matrix_isomorphism_compose_TT(self, x, y)
            X = self.T.toMatrix(x);
            Y = self.T.toMatrix(y);
            XY = X*Y;
            XY1 = self.T.toMatrix(self.T.compose(x, y));
            self.assert(isequal(XY, XY1));
        end
        function law_matrix_isomorphism_T(self, x)
            X = self.T.toMatrix(x);
            x1 = self.T.fromMatrix(X);
            self.T.assertEqv(x, x1);
        end
    end
end
