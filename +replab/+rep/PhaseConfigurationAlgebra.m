classdef PhaseConfigurationAlgebra < replab.rep.Algebra
% Algebra of matrices invariant under a monomial representation
    properties (SetAccess = protected)
        description;
        phaseConfiguration;
    end
    
    methods
        
        function self = PhaseConfigurationAlgebra(phaseConfiguration, fibers)
            self.n = phaseConfiguration.n;
            self.phaseConfiguration = phaseConfiguration;
            self.fibers = fibers;
            self.description = sprintf('Algebra of %d x %d matrices of dimension %d', self.n, self.n, self.phaseConfiguration.nOrbits);
        end
        
        function A = restricted(self, fiber)
            r = self.fibers.blockIndex(fiber);
            u = unique(r);
            newBlockIndex = zeros(1, length(fiber));
            for i = 1:length(u)
                j = u(i);
                newBlockIndex(r == j) = i;
            end
            newFibers = replab.Partition.fromBlockIndices(newBlockIndex);
            newPC = self.phaseConfiguration.restrict(fiber);
            A = replab.rep.PhaseConfigurationAlgebra(newPC, newFibers);
        end
        
        function M = sample(self)
            M = self.phaseConfiguration.sampleRealGaussian;
        end
        
        function M = sampleSelfAdjoint(self)
            M = self.phaseConfiguration.sampleSymmetricGaussian;
        end
        
        function T = project(self, T)
            T = self.phaseConfiguration.project(T);
        end
        
    end
    
end
