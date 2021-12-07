classdef ActionLaws < replab.Laws

    properties (SetAccess = protected)
        action % (`+replab.Action`): Action
        G % (`+replab.Group`): Group acting
        P % (`+replab.Domain`): Set acted upon
    end

    methods

        function self = ActionLaws(action)
            self.action = action;
            self.G = action.G;
            self.P = action.P;
        end

        function law_leftAction_compose_GGP(self, g, h, p)
        % Left action composition law
            gh = self.G.compose(g, h);
            g_h_p = self.action.leftAction(g, self.action.leftAction(h, p));
            gh_p = self.action.leftAction(gh, p);
            self.P.assertEqv(g_h_p, gh_p);
        end

        function law_rightAction_compose_PGG(self, p, g, h)
        % Right action composition law
            gh = self.G.compose(g, h);
            p_g_h = self.action.rightAction(self.action.rightAction(p, g), h);
            p_gh = self.action.rightAction(p, gh);
            self.P.assertEqv(p_g_h, p_gh);
        end

        function law_identity_leftAction_P(self, p)
        % Left action identity law
            self.P.assertEqv(p, self.action.leftAction(self.G.identity, p));
        end

        function law_identity_rightAction_P(self, p)
        % Right action identity law
            self.P.assertEqv(p, self.action.rightAction(p, self.G.identity));
        end

        function law_leftAction_rightAction_compatibility_GP(self, g, p)
        % Left/right action compatibility law
            gInv = self.G.inverse(g);
            self.P.assertEqv(self.action.leftAction(g, p), self.action.rightAction(p, gInv));
        end

    end

end
