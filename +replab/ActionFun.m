classdef ActionFun < replab.Action & replab.StrFun
    properties (SetAccess = protected)
        leftActionFun;
    end
    methods
        
        function self = ActionFun(description, G, P, leftActionFun)
            self = self@replab.StrFun(@(s, mc) description);
            self.G = G;
            self.P = P;
            self.leftActionFun = leftActionFun;
        end
    
    end
    
end
