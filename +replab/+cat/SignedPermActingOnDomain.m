classdef SignedPermActingOnDomain < replab.cat.Action
    
    methods
        
        function self = SignedPermActingOnDomain(n)
            self.groupCat = replab.cat.SignedPermAsGroup(n);
        end
        
        function i1 = leftAction(self, g, i)
            i1 = g(abs(i)) * sign(i);
        end
        
        function i1 = rightAction(self, i, g)
            i1 = find(abs(g) == abs(i));
            i1 = i1 * sign(g(i1)) * sign(i);
        end
        
    end
    
end