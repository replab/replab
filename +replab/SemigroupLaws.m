classdef SemigroupLaws < replab.DomainLaws
    methods
        function self = SemigroupLaws(T)
            self = self@replab.DomainLaws(T);
        end
        function law_associativity_TTT(self, x, y, z)
        % Checks associativity of group binary operation
            xy = self.T.compose(x, y);
            yz = self.T.compose(y, z);
            self.T.assertEqv(self.T.compose(xy, z), self.T.compose(x, yz));
        end
        function law_composeN_positive_TN10(self, x, n)
            xn1 = x;
            for i = 2:n
                xn1 = self.T.compose(xn1, x);
            end
            if n > 1
                xn2 = self.T.compose(self.T.composeN(x, n - 1), x);
                self.T.assertEqv(xn1, xn2);
            end
            xn3 = self.T.composeN(x, n);
            self.T.assertEqv(xn1, xn3);
        end
        function law_composeAll_TTT(self, x, y, z)
            xyz1 = self.T.compose(self.T.compose(x, y), z);
            xyz2 = self.T.composeAll({x y z});
            self.T.assertEqv(xyz1, xyz2);
        end
    end
end
