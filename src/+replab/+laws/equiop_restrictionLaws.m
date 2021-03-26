classdef equiop_restrictionLaws < replab.Laws

    properties (SetAccess = protected)
        G % (`+replab.CompactGroup`): Target group
    end

    methods

        function self = equiop_restrictionLaws(op)
            self@replab.laws.equiopLaws(op)
            self.G = op.target.group;
        end

        function law_representation_is_restriction_G(self, g)
            [t1, d1] = self.T.repR.image(g);
            [t2, d2] = self.T.repC.image(g);
            h = self.op.mu.imageElement(g);
            [s1, e1] = self.S.repR.image(g);
            [s2, e2] = self.S.repC.image(g);
            self.assertApproxEqual(s1, t1, d1+e1);
            self.assertApproxEqual(s2, t2, d2+e2);
        end

    end

end
