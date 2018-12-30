classdef RealConfigurationAlgebra < replab.RealCentralizerAlgebra
    
    properties (SetAccess = protected)
        phaseConfiguration;
    end
    
    methods

        function self = RealConfigurationAlgebra(realRep, phaseConfiguration)
            self = self@replab.RealCentralizerAlgebra(realRep);
            self.phaseConfiguration = phaseConfiguration;
        end
        
        function M = sampleGeneric(self)
            M = self.phaseConfiguration.sampleRealGaussian;
        end
        
        function M = sampleGenericSelfAdjoint(self)
            M = self.phaseConfiguration.sampleSymmetricGaussian;
        end
        
        function M1 = project(self, M)
            M1 = self.phaseConfiguration.project(M);
        end
        
    end
    
    methods (Static)
       
        function A = fromRealRep(realRep)
            nG = realRep.group.nGenerators;
            d = realRep.dimension;
            signedPerms = zeros(nG, d);
            for i = 1:nG
                signedPerms(i,:) = replab.SignedPermutations.fromMatrix(realRep.images{i});
            end
            pc = replab.rep.SignedConfigurationBuilder(signedPerms).toPhaseConfiguration;
            A = replab.RealConfigurationAlgebra(realRep, pc);
        end
        
    end
        
end
