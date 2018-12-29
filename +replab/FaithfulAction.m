classdef FaithfulAction < replab.Action
    
    methods % Abstract

        function p = findMovedElement(self, g)
        % Returns either p such that leftAction(g, p) != p
        % or [] if no such p exists
            f = self.findMovedElementFun;
            p = f(g);
        end
        
    end
    
end
