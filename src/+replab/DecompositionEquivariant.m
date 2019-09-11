classdef DecompositionEquivariant < replab.Equivariant
    
    methods

        function self = DecompositionEquivariant(repR, repC)
            replab.Dispatch.assert(isa(repR.group, 'replab.FiniteGroup'));
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
        % Projects any nR x nC matrix in the equivariant subspace
            assert(isa(self.group, 'replab.FiniteGroup'));
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                X = self.averageOver(X, T{i});
            end
        end
        
    end
    
end
