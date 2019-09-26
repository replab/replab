classdef CommutantFromEquivariant < replab.Commutant
% Derives the commutant from the equivariant space
    
    properties (Access = protected)
        equivariant_ % replab.Equivariant: Underlying equivariant space
    end
    
    methods (Access = protected)
        
        function e = equivariant(self)            
        % Equivariant space isomorphic to the commutant algebra as a vector space
        %
        % To avoid code duplication, we delegate computations involving the commutant to
        % the equivariant space of a representation to itself.
            if isempty(self.equivariant_)
                self.equivariant_ = replab.makeEquivariant(self.rep, self.rep);
            end
            e = self.equivariant_;
        end
        
    end

    methods
        
        function X = project(self, X)
        % Projects any n x n matrix in the invariant subspace
            X = self.equivariant.project(X);
        end

    end
    
end
