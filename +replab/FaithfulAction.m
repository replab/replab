classdef FaithfulAction < replab.Action
    
    methods % Abstract

        function p = findMovedElement(self, g)
        % Returns either p such that leftAction(g, p) != p
        % or [] if no such p exists
            error('Not implemented');
        end
        
    end

    methods (Static)
        
        function action = lambda(header, G, P, ...
                                 leftActionFun, findMovedElementFun)
            action = replab.lambda.FaithfulAction(header, G, P, ...
                                                  leftActionFun, findMovedElementFun);
        end
        
    end
    
end
