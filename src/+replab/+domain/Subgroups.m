classdef Subgroups < replab.domain.BoundedLattice
% Represents the lattice of subgroups of a finite group
%
% This class is mainly used to test the `+replab.FiniteGroup.closure` and `+replab.FiniteGroup.intersection` methods

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
            G = self.one.randomSubgroup;
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
