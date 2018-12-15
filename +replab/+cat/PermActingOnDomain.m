classdef PermActingOnDomain < replab.cat.Action
    
    methods
        
        function self = PermActingOnDomain(n)
            self.groupCat = replab.cat.PermAsGroup(n);
        end
        
        function i1 = leftAction(self, g, i)
            i1 = g(i);
        end
        
        function i1 = rightAction(self, i, g)
            i1 = find(g == i);
        end
        
    end
    
end