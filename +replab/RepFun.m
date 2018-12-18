classdef RepFun < replab.Rep & replab.cat.GroupMorphismFun
    
    methods
        
        function self = RepFun(group, fun, S, T)
            self@replab.cat.GroupMorphismFun(fun, S, T);
            self.group = group;
        end
                        
    end
    
end
