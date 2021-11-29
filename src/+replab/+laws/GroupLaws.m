classdef GroupLaws < replab.laws.MonoidLaws
% Group axioms

    methods

        function self = GroupLaws(S)
            self@replab.laws.MonoidLaws(S);
        end

        function law_composeN_integers_SN10(self, x, n)
        % Checks negative powers of a group element
            pow1 = self.S.composeN(x, n);
            pow2 = self.S.inverse(self.S.composeN(x, -n));
            self.S.assertEqv(pow1, pow2);
        end

        function law_inverse_S(self, x)
        % Checks that the composition with inverse is the identity
            xI = self.S.inverse(x);
            id1 = self.S.compose(x, xI);
            id2 = self.S.compose(xI, x);
            self.assert(self.S.isIdentity(id1));
            self.assert(self.S.isIdentity(id2));
        end

        function law_composeWithInverse_SS(self, x, y)
        % Checks that the composition with inverse method is properly implemented
            xyI = self.S.compose(x, self.S.inverse(y));
            self.S.assertEqv(xyI, self.S.composeWithInverse(x, y));
        end

        function law_leftConjugate_SS(self, x, y);
        % Checks that the left conjugate is correctly implemented
            xyxI1 = self.S.compose(x, self.S.composeWithInverse(y, x));
            xyxI2 = self.S.leftConjugate(x, y);
            self.S.assertEqv(xyxI1, xyxI2);
        end

        function law_inverse_compatible_with_compose_SS(self, x, y)
        % Checks that the inverse of a composition is the composition of swapped inverses
            xy = self.S.compose(x, y);
            yIxI = self.S.compose(self.S.inverse(y), self.S.inverse(x));
            self.S.assertEqv(self.S.inverse(xy), yIxI);
        end

    end

end
