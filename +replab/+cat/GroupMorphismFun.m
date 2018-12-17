classdef GroupMorphismFun < replab.cat.GroupMorphism
    
    properties
        fun; % morphism
    end
    
    methods
        
        function self = GroupMorphismFun(fun, S, T)
            self.fun = fun;
            self.S = S;
            self.T = T;
        end
        
        function t = image(self, s)
            t = self.fun(s);
        end
    
    end
    
end

