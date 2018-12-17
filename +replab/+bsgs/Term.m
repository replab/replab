classdef Term < replab.bsgs.Chain
% The terminal node in a BSGS chain
    
   methods
       
       function self = Term(A)
           self.A = A;
       end
       
       function s = sift(self, remaining)
           s = remaining;
       end
       
       function b = isTerm(self)
           b = true;
       end
       
   end
   
end
