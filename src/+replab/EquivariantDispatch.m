classdef EquivariantDispatch < replab.Dispatch
    
    properties (Constant)
        instance = replab.EquivariantDispatch;
    end
    
    methods
        
        function self = EquivariantDispatch
            self.register('FiniteGroup', @(repR, repC) replab.DecompositionEquivariant(repR, repC), 0);
        end
        
    end
    
end
