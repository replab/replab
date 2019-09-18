classdef SubRepNULaws < replab.RepLaws
    
    methods
        
        function self = SubRepNULaws(rep)
            self = self@replab.RepLaws(rep);
        end
        
        function law_basis_is_consistent(self)
            self.M.assertEqv(self.rep.F * self.rep.H, eye(self.rep.dimension));
        end
        
        function law_projector(self)
            proj = self.rep.projector;
            self.assert(~replab.isNonZeroMatrix(proj*proj - proj, replab.Settings.doubleEigTol));
        end
        
        function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                proj = self.rep.projector;
                rho = self.rep.image(g);
                self.assert(~replab.isNonZeroMatrix(proj*parentRho - parentRho*proj, replab.Settings.doubleEigTol));
                self.M.assertEqv(self.rep.F*parentRho*self.rep.H, rho);
            end
        end
        
    end
    
end
