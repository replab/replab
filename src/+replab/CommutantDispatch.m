classdef CommutantDispatch < replab.Dispatch
    
    properties (Constant)
        instance = replab.CommutantDispatch;
    end
    
    methods
        
        function self = CommuntantDispatch
            self.register('DecompositionCommutant', @(r) replab.rep.DecompositionCommutant(r), 5);
            self.register('ApproximateCommutant', @(r) replab.rep.ApproximateCommutant(r), 0);
        end
        
    end
    
end
