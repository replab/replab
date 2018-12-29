classdef Action < replab.Str
% A group action describing the action of elements of type G upon elements of type P
    properties (SetAccess = protected)
        G; % group
        P; % set acted upon
    end
    
    methods % Abstract methods
        
        function p1 = leftAction(self, g, p)
        % Returns the left action p1 = g(p) of G over P, which
        % is compatible with group composition in this way
        % p2 = g(h(p)) implies p2 = (g compose h)(p)
            f = self.leftActionFun;
            p1 = f(g, p);
        end
        
    end
    
    methods % Methods with default implementations
        
        function p1 = rightAction(self, p, g)
        % Returns the right action p1 = p^g, compatible with the
        % group composition as p2 = (p^g)^h = p^(g compose h)
            p1 = self.leftAction(self.G.inverse(g), p);
        end

    end
end
