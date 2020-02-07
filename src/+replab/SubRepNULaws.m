classdef SubRepNULaws < replab.RepLaws

    methods

        function self = SubRepNULaws(rep)
            self = self@replab.RepLaws(rep);
        end

        function law_image_relation_with_parent_rep_G_M(self, g, m)
            m1 = self.rep.H * self.rep.image(g) * m;
            m2 = self.rep.parent.image(g) * self.rep.H * m;
            self.assert(~replab.isNonZeroMatrix(m1 - m2, replab.Parameters.doubleEigTol));
        end

    end

end
