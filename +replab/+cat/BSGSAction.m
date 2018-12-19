classdef BSGSAction < replab.cat.Action
    
    methods

        function p = findMovedElement(self, g)
        % Returns either p such that leftAction(g, p) != p
        % or [] if no such p exists
            f = self.findMovedElementFun;
            p = f(g);
        end
        
    end
    
    methods % LAWS
        
        function law_findMovedElement_identity(self)
            p = self.findMovedElement(self.G.identity);
            assertEqual(p, []);
        end
        
        function law_findMovedElement_G(self, g)
            p = self.findMovedElement(g);
            if ~isequal(p, [])
                p1 = self.leftAction(g, p);
                self.P.assertNotEqv(p, p1);
            end
        end
        
    end
    
end
