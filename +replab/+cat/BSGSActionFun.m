classdef BSGSActionFun < replab.cat.ActionFun & replab.cat.BSGSAction
    
    properties (SetAccess = protected)
        findMovedElementFun;
    end
    
    methods
        
        function self = BSGSActionFun(description, G, P, leftActionFun, findMovedElementFun)
            self@replab.cat.ActionFun(description, G, P, leftActionFun);
            self.findMovedElementFun = findMovedElementFun;
        end
        
        function p = findMovedElement(self, g)
            p = self.findMovedElementFun(g);
        end
        
    end
    
end
