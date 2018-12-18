classdef ActionFun < replab.cat.Action
% 
    properties (SetAccess = protected)
        description;
        leftActionFun;
    end
    
    methods
        
        function self = ActionFun(description, G, P, leftActionFun)
            self.description = description;
            self.G = G;
            self.P = P;
            self.leftActionFun = leftActionFun;
        end
        
        function s = str(self)
            s = self.description;
        end
        
        function p1 = leftAction(self, g, p)
            p1 = self.leftActionFun(g, p);
        end

    end    
    
end
