classdef equiop_restrictionLaws < replab.laws.equiopLaws

    properties (SetAccess = protected)
        G % (`+replab.CompactGroup`): Target group
    end

    methods

        function self = equiop_restrictionLaws(op)
            self@replab.laws.equiopLaws(op)
            self.G = op.target.group;
        end

        function law_representation_is_restriction_G(self, g)
            d1 = self.T.repR.errorBound;
            d2 = self.T.repC.errorBound;
            t1 = self.T.repR.image(g);
            t2 = self.T.repC.image(g);
            h = self.op.mu.imageElement(g);
            e1 = self.S.repR.errorBound;
            e2 = self.S.repC.errorBound;
            s1 = self.S.repR.image(h);
            s2 = self.S.repC.image(h);
            self.assertApproxEqual(s1, t1, d1+e1);
            self.assertApproxEqual(s2, t2, d2+e2);
        end

    end

end
