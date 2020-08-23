classdef SubRepLaws < replab.laws.RepLaws
% Laws for subrepresentations

    methods

        function self = SubRepLaws(rep)
            self = self@replab.laws.RepLaws(rep);
        end

        function law_basis_and_internal_embedding_(self)
            self.M.assertEqv(self.rep.E_internal * self.rep.B_internal, eye(self.rep.dimension));
            P = self.rep.B_internal * self.rep.E_internal;
            self.M.assertEqv(P*P, P);
        end

        function law_image_relation_with_parent_rep_GM(self, g, m)
            m1 = full(self.rep.B_internal * self.rep.image(g) * m);
            m2 = full(self.rep.parent.image(g) * self.rep.B_internal * m);
            self.assert(~replab.isNonZeroMatrix(m1 - m2, replab.Parameters.doubleEigTol));
        end

         function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                rho = self.rep.image(g);
                self.M.assertEqv(full(self.rep.E_internal*parentRho*self.rep.B_internal), rho);
            end
        end

    end

end
