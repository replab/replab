classdef SubRepLaws < replab.RepLaws
    
    methods
        
        function self = SubRepLaws(rep)
            self = self@replab.RepLaws(rep);
        end
        
        function law_basis_is_unitary(self)
            self.M.assertEqv(self.rep.U * self.rep.U', eye(self.rep.dimension));
        end
        
        function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                proj = self.rep.projector;
                rho = self.rep.image(g);
                self.assert(~replab.isNonZeroMatrix(proj*parentRho - parentRho*proj, replab.Settings.doubleEigTol));
                self.M.assertEqv(self.rep.U*parentRho*self.rep.U', rho);
            end
        end
        
        function law_nice_basis_reproduces_basis(self)
            if ~isempty(self.rep.niceBasis)
                self.assert(~replab.isNonZeroMatrix(self.rep.niceBasis.U - self.rep.U, replab.Settings.doubleEigTol));
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
