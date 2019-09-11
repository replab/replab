classdef FaithfulActionLaws < replab.ActionLaws
    methods
        function self = FaithfulActionLaws(action)
            self@replab.ActionLaws(action);
        end
        function law_findMovedElement_identity(self)
            p = self.action.findMovedElement(self.G.identity);
            assertEqual(p, []);
        end
        function law_findMovedElement_G(self, g)
            p = self.action.findMovedElement(g);
            if ~isempty(p)
                p1 = self.action.leftAction(g, p);
                self.P.assertNotEqv(p, p1);
            end
        end
    end
end
