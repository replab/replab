classdef DecompositionDispatch < replab.Dispatch

    properties (Constant)
        instance = replab.rep.DecompositionDispatch;
    end
    
    methods
        
        function self = DecompositionDispatch
            self.register('Generic', @replab.rep.decomposition, 0);
            self.register('Table-based symmetric group', @replab.rep.decomposeSymmetric, 10)
        end
    end
    
end
