classdef DecompositionCommutant < replab.Commutant
% Commutant algebra of a finite group representation 
    
    properties (SetAccess = protected)
        
    end
    
    methods
        
        function self = DecompositionCommutant(rep)
            self@replab.Commutant(rep);
            assert(isa(rep.group, 'replab.FiniteGroup'));
        end
        
        function X = project(self, X)
        % Projects any n x n matrix in the equivariant subspace
            T = self.group.decomposition.transversals;
            for i = length(T):-1:1
                X = self.averageOver(X, T{i});
            end
        end
        
    end
    
end
