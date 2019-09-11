classdef ApproximateCommutant < replab.Commutant
% Commutant algebra of a finite group representation 
    
    methods
        
        function self = ApproximateCommutant(rep)
            self@replab.Commutant(rep);
        end
        
        function X = project(self, X)
        % Projects any n x n matrix in the equivariant subspace
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                X = self.averageOver(X, T{i});
            end
        end
        
    end
    
end
