classdef PhaseConfigurationAlgebra < replab.rep.Algebra

    properties (SetAccess = protected)
        description;
        phaseConfiguration;
    end
    
    methods
        
        function self = PhaseConfigurationAlgebra(phaseConfiguration, partition)
            self.n = phaseConfiguration.n;
            self.phaseConfiguration = phaseConfiguration;
            self.description = sprintf('Algebra of %d x %d matrices of dimension %d', self.n, self.n, self.phaseConfiguration.nOrbits);
            self.partition = partition;
        end
        
        function M = sample(self)
            M = self.phaseConfiguration.sampleRealGaussian;
        end
        
        function M = sampleSelfAdjoint(self)
            M = self.phaseConfiguration.sampleSymmetricGaussian;
        end
        
        function M = project(self, T)
            T = self.phaseConfiguration.project(T);
        end
        
    end
    
end
