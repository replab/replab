classdef Subgroups < replab.domain.BoundedLattice

    methods

        function self = Subgroups(group)
            self@replab.domain.BoundedLattice(group.trivialSubgroup, group);
        end

    end

    methods % Implementations

        function b = eqv(self, lhs, rhs)
            b = (lhs == rhs);
        end

        function G = sample(self)
            try
                G = self.one.randomProperSubgroup;
            catch
                if randi(2) == 1
                    G = self.one;
                else
                    G = self.zero;
                end
            end
        end

        function b = ltEqv(self, lhs, rhs)
            b = lhs.isSubgroupOf(rhs);
        end

        function z = join(self, x, y)
            z = x.closure(y);
        end

        function z = meet(self, x, y)
            z = x.intersection(y);
        end

    end

end
