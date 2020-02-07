classdef SubRepLaws < replab.RepLaws

    methods

        function self = SubRepLaws(rep)
            self = self@replab.RepLaws(rep);
        end

        function law_basis_is_unitary_(self)
            self.M.assertEqv(self.rep.U * self.rep.U', eye(self.rep.dimension));
        end

        function law_image_relation_with_parent_rep_GM(self, g, m)
            m1 = self.rep.U' * self.rep.image(g) * m;
            m2 = self.rep.parent.image(g) * self.rep.U' * m;
            self.assert(~replab.isNonZeroMatrix(m1 - m2, replab.Parameters.doubleEigTol));
        end

        function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                rho = self.rep.image(g);
                self.M.assertEqv(self.rep.U*parentRho*self.rep.U', rho);
            end
        end

        function law_nice_basis_reproduces_basis_(self)
            if ~isempty(self.rep.niceBasis)
                self.assert(~replab.isNonZeroMatrix(self.rep.niceBasis.U - self.rep.U, replab.Parameters.doubleEigTol));
            end
        end

        function law_respects_division_algebra_G(self, g)
            if isequal(self.rep.field, 'R') && self.rep.isKnownCanonicalIrreducible
                rho = self.rep.image(g);
                if isequal(self.rep.irrepInfo.divisionAlgebra, 'C')
                    rho1 = replab.domain.ComplexTypeMatrices.project(rho);
                    self.M.assertEqv(rho, rho1);
                elseif isequal(self.rep.irrepInfo.divisionAlgebra, 'H')
                    rho1 = replab.domain.QuaternionTypeMatrices.project(rho);
                    self.M.assertEqv(rho, rho1);
                end
            end
        end

    end

end
