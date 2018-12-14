classdef Rep < handle
    
    properties (SetAccess = protected)
        permGrp;
        d; % Representation dimension
        field; % either 'R' or 'C'
    end
    
    methods
       
        function M = image(self, permutation)
            error('Not implemented');
        end
        
        %rho = sampleGroupElement
        
        %rho = sampleGroupAlgebra

        %rho1 = projectInvariantMatrix(rho)
        
        %rho1 = projectSelfAdjointInvariantMatrix(rho)
        
        
        %P = coordinateOrbits
        
    end
    
end
