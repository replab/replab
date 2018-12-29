classdef FaithfulActionFun < replab.ActionFun & replab.FaithfulAction
    
    properties (SetAccess = protected)
        findMovedElementFun;
    end
    
    methods
        
        function self = FaithfulActionFun(description, G, P, leftActionFun, findMovedElementFun)
            self@replab.ActionFun(description, G, P, leftActionFun);
            self.findMovedElementFun = findMovedElementFun;
        end
        
    end
    
end
