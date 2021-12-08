classdef BoundedLatticeLaws < replab.laws.DomainLaws

    methods

        function self = BoundedLatticeLaws(S)
            self@replab.laws.DomainLaws(S);
        end

        function law_meet_associativity_SSS(self, x, y, z)
            xy = self.S.meet(x, y);
            yz = self.S.meet(y, z);
            self.S.assertEqv(self.S.meet(xy, z), self.S.meet(x, yz));
        end

        function law_meet_commutative_SS(self, x, y)
            self.S.assertEqv(self.S.meet(x, y), self.S.meet(y, x));
        end

        function law_meet_idempotent_S(self, x)
            self.S.assertEqv(self.S.meet(x, x), x);
        end

        function law_join_associativity_SSS(self, x, y, z)
            xy = self.S.join(x, y);
            yz = self.S.join(y, z);
            self.S.assertEqv(self.S.join(xy, z), self.S.join(x, yz));
        end

        function law_join_commutative_SS(self, x, y)
            self.S.assertEqv(self.S.join(x, y), self.S.join(y, x));
        end

        function law_join_idempotent_S(self, x)
            self.S.assertEqv(self.S.join(x, x), x);
        end

        function law_absorption_SS(self, x, y)
            self.S.assertEqv(self.S.join(x, self.S.meet(x, y)), self.S.meet(x, self.S.join(x, y)));
        end

        function law_bounded_join_S(self, x)
            self.S.assertEqv(self.S.join(self.S.zero, x), x);
            self.S.assertEqv(self.S.join(x, self.S.zero), x);
        end

        function law_bounded_meet_S(self, x)
            self.S.assertEqv(self.S.meet(self.S.one, x), x);
            self.S.assertEqv(self.S.meet(x, self.S.one), x);
        end

        function law_partial_order_meet_SS(self, x, y)
            assert(self.S.ltEqv(x, y) == self.S.eqv(self.S.meet(x, y), x));
        end

        function law_partial_order_join_SS(self, x, y)
            assert(self.S.ltEqv(x, y) == self.S.eqv(self.S.join(x, y), y));
        end

    end

end
