classdef Algebra < replab.Str
% A n x n real matrix algebra
    
    properties (SetAccess = protected)
        fibers; % A fiber is a subset of the coordinates 1..n
                % such that { M(fiber, fiber) for M \in algebra }
                % is also an algebra
                % The finest decomposition of 1..n into fibers
                % is stored as a replab.Partition into this variable
        n; % Size of the n x n matrix representation of this algebra
    end
    
    methods
        
        function A = restricted(self, fiber)
        % Returns a restriction of this algebra to the given fiber
            error('Not implemented');
        end
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
            phaseConfiguration = replab.rep.PhaseConfiguration.fromSignedPerm(signedPerms);
            A = replab.rep.PhaseConfigurationAlgebra(phaseConfiguration, rep.fibers);
        end
        
    end
    
end
