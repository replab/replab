classdef RepFun < replab.Rep & replab.cat.GroupMorphismFun
    
    methods
        
        function self = RepFun(fun, S, T)
            self@replab.cat.GroupMorphismFun(fun, S, T);
        end
                        
    end
    
end
