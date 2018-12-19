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
        
    end    
    
end
