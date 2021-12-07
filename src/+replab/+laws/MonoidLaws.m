classdef MonoidLaws < replab.laws.DomainLaws
% Laws obeyed by a monoid

    properties (SetAccess = protected)
        N10 % (`+replab.Domain`): Integer between 1 and 10
    end

    methods

        function self = MonoidLaws(S)
            self@replab.laws.DomainLaws(S);
            self.N10 = replab.domain.intAsDouble(1, 10);
        end

        function law_associativity_SSS(self, x, y, z)
        % Checks associativity of group binary operation
            xy = self.S.compose(x, y);
            yz = self.S.compose(y, z);
            self.S.assertEqv(self.S.compose(xy, z), self.S.compose(x, yz));
        end

        function law_composeN_positive_SN10(self, x, n)
        % Checks repeated composition of an element with itself
            xn1 = x;
            for i = 2:n
                xn1 = self.S.compose(xn1, x);
            end
            if n > 1
                xn2 = self.S.compose(self.S.composeN(x, n - 1), x);
                self.S.assertEqv(xn1, xn2);
            end
            xn3 = self.S.composeN(x, n);
            self.S.assertEqv(xn1, xn3);
        end

        function law_composeAll_SSS(self, x, y, z)
        % Checks composition of a list of items in a cell array
            xyz1 = self.S.compose(self.S.compose(x, y), z);
            xyz2 = self.S.composeAll({x y z});
            self.S.assertEqv(xyz1, xyz2);
        end

        function law_composeN_zero_S(self, x)
        % Checks composition of zero times an element with itself
            id1 = self.S.composeN(x, 0);
            id2 = self.S.identity;
            self.S.assertEqv(id1, id2);
        end

        function law_identity_(self)
        % Checks that the identity element is equal to the identity
            id = self.S.identity;
            self.assert(self.S.isIdentity(id));
        end

    end

end
