classdef SubRepLaws < replab.RepLaws

    methods

        function self = SubRepLaws(rep)
            self = self@replab.RepLaws(rep);
        end

        function law_basis_and_internal_embedding_(self)
            self.M.assertEqv(self.rep.F_internal * self.rep.H_internal, eye(self.rep.dimension));
        end

        function law_image_relation_with_parent_rep_GM(self, g, m)
            m1 = full(self.rep.H_internal * self.rep.image(g) * m);
            m2 = full(self.rep.parent.image(g) * self.rep.H_internal * m);
            self.assert(~replab.isNonZeroMatrix(m1 - m2, replab.Parameters.doubleEigTol));
        end

         function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                rho = self.rep.image(g);
                self.M.assertEqv(full(self.rep.F_internal*parentRho*self.rep.H_internal), rho);
            end
        end

    end

end
