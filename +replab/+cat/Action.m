classdef Action < replab.cat.Laws
    
    properties (SetAccess = protected)
        G; % group
        P; % set acted upon
    end
    
    methods
        
        % Implement
        %
        % function p1 = leftAction(self, g, p)

        function p1 = rightAction(self, p, g)
            p1 = self.leftAction(self.G.inverse(g), p);
        end

    end
    
    methods % LAWS
        
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
