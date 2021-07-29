classdef BoundedLatticeLaws < replab.laws.DomainLaws

    methods

        function self = BoundedLatticeLaws(T)
            self@replab.laws.DomainLaws(T);
        end

        function law_meet_associativity_TTT(self, x, y, z)
            xy = self.T.meet(x, y);
            yz = self.T.meet(y, z);
            self.T.assertEqv(self.T.meet(xy, z), self.T.meet(x, yz));
        end

        function law_meet_commutative_TT(self, x, y)
            self.T.assertEqv(self.T.meet(x, y), self.T.meet(y, x));
        end

        function law_meet_idempotent_T(self, x)
            self.T.assertEqv(self.T.meet(x, x), x);
        end

        function law_join_associativity_TTT(self, x, y, z)
            xy = self.T.join(x, y);
            yz = self.T.join(y, z);
            self.T.assertEqv(self.T.join(xy, z), self.T.join(x, yz));
        end

        function law_join_commutative_TT(self, x, y)
            self.T.assertEqv(self.T.join(x, y), self.T.join(y, x));
        end

        function law_join_idempotent_T(self, x)
            self.T.assertEqv(self.T.join(x, x), x);
        end

        function law_absorption_TT(self, x, y)
            self.T.assertEqv(self.T.join(x, self.T.meet(x, y)), self.T.meet(x, self.T.join(x, y)));
        end

        function law_bounded_join_T(self, x)
            self.T.assertEqv(self.T.join(self.T.zero, x), x);
            self.T.assertEqv(self.T.join(x, self.T.zero), x);
        end

        function law_bounded_meet_T(self, x)
            self.T.assertEqv(self.T.meet(self.T.one, x), x);
            self.T.assertEqv(self.T.meet(x, self.T.one), x);
        end

        function law_partial_order_meet_TT(self, x, y)
            assert(self.T.ltEqv(x, y) == self.T.eqv(self.T.meet(x, y), x));
        end

        function law_partial_order_join_TT(self, x, y)
            assert(self.T.ltEqv(x, y) == self.T.eqv(self.T.join(x, y), y));
        end

    end

end
