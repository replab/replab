classdef DecompositionDispatch < replab.Dispatch

    properties (Constant)
        instance = replab.rep.DecompositionDispatch;
    end
    
    methods
        
        function self = DecompositionDispatch
            self.register('Generic', @replab.rep.decomposition, 0);
        end
        
    end
    
end
