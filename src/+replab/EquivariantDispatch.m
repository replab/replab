classdef EquivariantDispatch < replab.Dispatch
    
    properties (Constant)
        instance = replab.EquivariantDispatch;
    end
    
    methods
        
        function self = EquivariantDispatch
            self.register('ForRepByImages', @(repR, repC) replab.equivariant.ForRepByImages(repR, repC), 10);
            self.register('ForFiniteGroup', @(repR, repC) replab.equivariant.ForFiniteGroup(repR, repC), 5);
            % Default method, works for all compact groups
            self.register('ForCompactGroup', @(repR, repC) replab.equivariant.ForCompactGroup(repR, repC), 0);
        end
        
    end
    
end
