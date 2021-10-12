classdef TotalOrderLaws < replab.laws.DomainLaws

    methods

        function self = TotalOrderLaws(T)
            self@replab.laws.DomainLaws(T);
        end

        function law_reflexive_T(self, x)
            self.assert(self.T.compare(x, x) == 0);
        end

        function law_transitive_TTT(self, x, y, z)
            if self.T.compare(x, y) <= 0 && self.T.compare(y, z) <= 0
                self.assert(self.T.compare(x, z) <= 0);
            end
        end

    end

end
