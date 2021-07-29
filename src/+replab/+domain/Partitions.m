classdef Partitions < replab.domain.BoundedLattice

    properties (SetAccess = protected)
        n % (integer): Underlying set size
    end

    methods

        function self = Partitions(n)
            self@replab.domain.BoundedLattice(replab.Partition.finest(n), replab.Partition.trivial(n));
            self.n = n;
        end

    end

    methods % Implementation

        function P = sample(self)
            if self.n == 1 || (self.n == 2 && randi(2) == 1)
                P = self.zero;
            elseif self.n == 2
                P = self.one;
            else
                switch randi(8)
                  case 1
                    P = self.zero;
                  case 2
                    P = self.one;
                  otherwise
                    P = self.zero;
                    while P.nBlocks == 1 || P.nBlocks == self.n
                        P = replab.Partition.fromVector(randi(3, 1, self.n));
                    end
                end
            end
        end

        function b = eqv(self, lhs, rhs)
            b = lhs == rhs;
        end

        function b = ltEqv(self, lhs, rhs)
            b = lhs <= rhs;
        end

        function z = join(self, x, y)
            z = x.join(y);
        end

        function z = meet(self, x, y)
            z = x.meet(y);
        end

    end

end
