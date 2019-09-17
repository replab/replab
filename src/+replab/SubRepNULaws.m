classdef SubRepNULaws < replab.nu.RepLaws
    
    methods
        
        function self = SubRepLaws(rep)
            self = self@replab.nu.RepLaws(rep);
        end
        
        function law_basis_is_consistent(self)
            self.M.assertEqv(self.rep.F * self.rep.H, eye(self.rep.dimension));
        end
        
        function law_projector(self)
            proj = self.rep.projector;
            self.assert(~replab.isNonZeroMatrix(proj*proj - proj));
        end
        
        function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                proj = self.rep.projector;
                rho = self.rep.image(g);
                self.assert(~replab.isNonZeroMatrix(proj*parentRho - parentRho*proj, replab.Parameters.doubleEigTol));
                self.M.assertEqv(self.rep.U*parentRho*self.rep.U', rho);
            end
        end
        
    end
    
end
