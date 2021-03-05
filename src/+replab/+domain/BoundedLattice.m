classdef BoundedLattice < replab.Domain
% Describes a bounded lattice
%
% See `<https://www.math24.net/lattices>_`
    properties (SetAccess = protected)
        zero % Greatest lower bound
        one % Least upper bound
    end

    methods

        function self = BoundedLattice(zero, one)
            self.zero = zero;
            self.one = one;
        end

    end

    methods % Implementation

        function l = laws(self)
            l = replab.laws.BoundedLatticeLaws(self);
        end

    end

    methods % Partial order operations

        function b = ltEqv(self, lhs, rhs)
            error('Abstract');
        end

        function b = gtEqv(self, lhs, rhs)
            b = self.ltEqv(rhs, lhs);
        end

    end

    methods % Lattice operations

        function z = join(self, x, y)
            error('Abstract');
        end

        function z = meet(self, x, y)
            error('Abstract');
        end

    end

end
