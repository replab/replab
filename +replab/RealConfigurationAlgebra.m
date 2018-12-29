classdef RealConfigurationAlgebra < replab.RealCentralizerAlgebra
    
    properties (SetAccess = protected)
        phaseConfiguration;
    end
    
    methods

        function self = RealConfigurationAlgebra(realRep, phaseConfiguration)
            self = self@replab.RealCentralizerAlgebra(realRep);
            self.phaseConfiguration = phaseConfiguration;
        end
        
        function M = sampleUniformly(self)
            M = self.phaseConfiguration.sampleRealGaussian;
        end
        
        function M = sampleUniformlySelfAdjoint(self)
            M = self.phaseConfiguration.sampleSymmetricGaussian;
        end
        
        function M1 = project(self, M)
            M1 = self.phaseConfiguration.project(M);
        end
    end
        
end
