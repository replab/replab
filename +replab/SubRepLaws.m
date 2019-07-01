classdef SubRepLaws < replab.RepLaws
    
    methods
        
        function law_relation_with_parent_rep_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                proj = self.rep.U*self.rep.U';
                rho = self.rep.image(g);
                self.assert(~replab.isNonZeroMatrix(proj*parentRho - parentRho*proj, replab.Settings.doubleEigTol));
                self.M.assertEqv(self.rep.U'*parentRho*self.rep.U, rho);
            end
        end
        
    end
    
end
