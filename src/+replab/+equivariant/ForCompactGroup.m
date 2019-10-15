classdef ForCompactGroup < replab.Equivariant

    methods

        function self = ForCompactGroup(repR, repC)
            self = self@replab.Equivariant(repR, repC);
        end

        function X = project(self, X)
            nI = replab.Parameters.averagingIterations;
            nS = replab.Parameters.averagingSamples;
            for i = 1:nI
                g = self.group.sample;
                S = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X));
                for j = 2:nS
                    g = self.group.sample;
                    S = S + self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X));
                end
                X = S/nS;
            end
        end

    end

end
