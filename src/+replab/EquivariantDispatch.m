classdef EquivariantDispatch < replab.Dispatch
    
    properties (Constant)
        instance = replab.EquivariantDispatch;
    end
    
    methods
        
        function self = EquivariantDispatch
            self.register('ForRepByImages', @(repR, repC) replab.equivariant.ForRepByImages(repR, repC), 10);
            self.register('ForFiniteGroup', @(repR, repC) replab.equivariant.ForFiniteGroup(repR, repC), 5);
        end
        
    end
    
end
