classdef GroupLaws < replab.MonoidLaws
% Group axioms

    methods

        function self = GroupLaws(T)
            self@replab.MonoidLaws(T);
        end

        function law_composeN_integers_TN10(self, x, n)
        % Checks negative powers of a group element
            pow1 = self.T.composeN(x, n);
            pow2 = self.T.inverse(self.T.composeN(x, -n));
            self.T.assertEqv(pow1, pow2);
        end

        function law_inverse_T(self, x)
        % Checks that the composition with inverse is the identity
            xI = self.T.inverse(x);
            id1 = self.T.compose(x, xI);
            id2 = self.T.compose(xI, x);
            self.assert(self.T.isIdentity(id1));
            self.assert(self.T.isIdentity(id2));
        end

        function law_composeWithInverse_TT(self, x, y)
        % Checks that the composition with inverse method is properly implemented
            xyI = self.T.compose(x, self.T.inverse(y));
            self.T.assertEqv(xyI, self.T.composeWithInverse(x, y));
        end

        function law_inverse_compatible_with_compose_TT(self, x, y)
        % Checks that the inverse of a composition is the composition of swapped inverses
            xy = self.T.compose(x, y);
            yIxI = self.T.compose(self.T.inverse(y), self.T.inverse(x));
            self.T.assertEqv(self.T.inverse(xy), yIxI);
        end

    end

end
