classdef FaithfulActionFun < replab.FaithfulAction
    
    properties (SetAccess = protected)
        header;
        leftActionFun;
        findMovedElementFun;
    end
    
    methods
        
        function self = FaithfulActionFun(header, G, P, leftActionFun, findMovedElementFun)
            self.header = header;
            self.G = G;
            self.P = P;
            self.leftActionFun = leftActionFun;
            self.findMovedElementFun = findMovedElementFun;
        end
        
        function str = headerStr(self)
            str = self.header;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.ActionFun(self), ...
                hiddenFields@replab.FaithfulAction(self), ...
                {'header'}, ...
                );
        end

        function p1 = leftAction(self, g, p)
            f = self.leftActionFun;
            p1 = f(g, p);
        end
        
        function p = findMovedElement(self, g)
            f = self.findMovedElementFun;
            p = f(g);
        end
        
    end
    
end
