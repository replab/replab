classdef Start < handle
    
    properties
        next;
    end
   
    methods (Static)
        
        function start = emptyChain(A)
            start = replab.bsgs.Start;
            start.next = replab.bsgs.Term(A);
        end
        
    end

end
