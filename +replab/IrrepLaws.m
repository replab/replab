classdef IrrepLaws < replab.SubRepLaws
    
    methods

        function self = IrrepLaws(rep)
            self = self@replab.SubRepLaws(rep);
        end
        
        function law_respects_division_algebra_G(self, g)
            if isequal(self.rep.field, 'R')
                rho = self.rep.image(g);
                rho1 = self.rep.realDivisionAlgebra.projectMatrix(rho);
                self.M.assertEqv(rho, rho1);
            end
        end
        
    end

end
