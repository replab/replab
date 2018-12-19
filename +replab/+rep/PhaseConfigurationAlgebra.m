classdef PhaseConfigurationAlgebra < replab.rep.Algebra
% Algebra of matrices invariant under a monomial representation
    properties (SetAccess = protected)
        description;
        phaseConfiguration;
    end
    
    methods
        
        function self = PhaseConfigurationAlgebra(phaseConfiguration)
            self.n = phaseConfiguration.n;
            self.phaseConfiguration = phaseConfiguration;
            self.description = sprintf('Algebra of %d x %d matrices of dimension %d', self.n, self.n, self.phaseConfiguration.nOrbits);
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
