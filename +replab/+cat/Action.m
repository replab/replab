classdef Action < replab.cat.Laws
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
    
    methods % Laws
        
        function law_leftAction_compose_GGP(self, g, h, p)
            gh = self.G.compose(g, h);
            g_h_p = self.leftAction(g, self.leftAction(h, p));
            gh_p = self.leftAction(gh, p);
            self.P.assertEqv(g_h_p, gh_p);
        end
        
        function law_rightAction_compose_PGG(self, p, g, h)
            gh = self.G.compose(g, h);
            p_g_h = self.rightAction(self.rightAction(p, g), h);
            p_gh = self.rightAction(p, gh);
            self.P.assertEqv(p_g_h, p_gh);
        end
        
        function law_identity_leftAction_P(self, p)
            self.P.assertEqv(p, self.leftAction(self.G.identity, p));
        end
        
        function law_identity_rightAction_P(self, p)
            self.P.assertEqv(p, self.rightAction(p, self.G.identity));
        end
        
        function law_leftAction_rightAction_compatibility_GP(self, g, p)
            gInv = self.G.inverse(g);
            self.P.assertEqv(self.leftAction(g, p), self.rightAction(p, gInv));
        end
        
    end
    
end
