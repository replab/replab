classdef MonoidLaws < replab.DomainLaws
% Laws obeyed by a monoid

    properties
        N10 % Integer between 1 and 10
    end

    methods

        function self = MonoidLaws(T)
            self@replab.DomainLaws(T);
            self.N10 = replab.domain.intAsDouble(1, 10);
        end

        function law_associativity_TTT(self, x, y, z)
        % Checks associativity of group binary operation
            xy = self.T.compose(x, y);
            yz = self.T.compose(y, z);
            self.T.assertEqv(self.T.compose(xy, z), self.T.compose(x, yz));
        end

        function law_composeN_positive_TN10(self, x, n)
        % Checks repeated composition of an element with itself
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
        % Checks composition of a list of items in a cell array
            xyz1 = self.T.compose(self.T.compose(x, y), z);
            xyz2 = self.T.composeAll({x y z});
            self.T.assertEqv(xyz1, xyz2);
        end

        function law_composeN_zero_T(self, x)
        % Checks composition of zero times an element with itself
            id1 = self.T.composeN(x, 0);
            id2 = self.T.identity;
            self.T.assertEqv(id1, id2);
        end

        function law_identity_(self)
        % Checks that the identity element is equal to the identity
            id = self.T.identity;
            self.assert(self.T.isIdentity(id));
        end

    end

end
