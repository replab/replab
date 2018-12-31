classdef ActionFun < replab.Action & replab.Str
    properties (SetAccess = protected)
        leftActionFun;
    end
    methods
        
        function self = ActionFun(description, G, P, leftActionFun)
            self = self@replab.Str(description);
            self.G = G;
            self.P = P;
            self.leftActionFun = leftActionFun;
        end
    end    
    
end
