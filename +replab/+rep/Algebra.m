classdef Algebra < replab.Str

    properties (SetAccess = protected)
        n;
        partition;
    end
    
    methods
        
        function M = sample(self)
            error('Not implemented');
        end
        
        function M = sampleSelfAdjoint(self)
            error('Not implemented');
        end
        
        function M = project(self, T)
            error('Not implemented');
        end
        
    end
    
    methods (Static)
       
        function A = forNonSignedPermRep(rep, matrices)
            error('TODO: not implemented.')
        end
        
        function A = forRep(rep)
            nG = rep.group.nGenerators;
            matrices = cell(1, nG);
            d = rep.dimension;
            for i = 1:nG
                matrices{i} = rep.image(rep.group.generator(i));
            end
            signedPerms = zeros(nG, d);
            for i = 1:nG
                sp = replab.SignedPermutations.fromMatrix(matrices{i});
                if isempty(sp)
                    A = replab.rep.Algebra.forNonSignedPermRep;
                    return
                end
                signedPerms(i,:) = sp;
            end
            partition = replab.Partition.permutationsOrbits(abs(signedPerms));
            phaseConfiguration = replab.rep.PhaseConfiguration.fromSignedPerm(signedPerms);
            A = replab.rep.PhaseConfigurationAlgebra(phaseConfiguration, partition);
        end
        
    end
    
end
