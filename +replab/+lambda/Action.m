classdef Action < replab.Action
    
    properties (SetAccess = protected)
        header;
        leftActionFun;
    end
    
    methods
        
        function self = Action(header, G, P, leftActionFun)
            self.header = header;
            self.G = G;
            self.P = P;
            self.leftActionFun = leftActionFun;
        end
        
        function str = headerStr(self)
            str = self.header;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Action(self), ...
                {'header'} ...
                );
        end
    
        function p1 = leftAction(self, g, p)
            f = self.leftActionFun;
            p1 = f(g, p);
        end
        
    end
    
end
