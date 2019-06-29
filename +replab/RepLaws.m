classdef RepLaws < replab.Laws
    properties
        rep;
        G; % group of which rep is a representation
           %C; % centralizer algebra
        M; % GL_d(R or C)
    end
    methods
        function self = RepLaws(rep)
            self.rep = rep;
            d = self.rep.dimension;
            self.G = rep.group;
            %self.C = rep.centralizerAlgebra;
            eqvFun = @(X, Y) ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
            switch rep.field
              case 'R'
                self.M = replab.GroupFun('Mat', eqvFun, @() replab.rep.sampleRealMatrix(d, d), ...
                                         @(x, y) x*y, eye(d), @(x) inv(x));
              case 'C'
                self.M = replab.GroupFun('Mat', eqvFun, @() replab.rep.sampleComplexMatrix(d, d), ...
                                         @(x, y) x*y, eye(d), @(x) inv(x));                
              otherwise
                error('Unknown field')
            end
            
        end
        function morphismLaws = laws_asGroupHomomorphism(self)
            morphismLaws = replab.GroupMorphismLaws(@(g) self.rep.image(g), self.G, self.R);
        end
        function law_commutes_with_centralizer_algebra_GC(self, g, C)
            rho = self.rep.image(g);
            self.R.assertEqv(rho*C, C*rho);
        end
        function law_parent_G(self, g)
            if ~isempty(self.rep.parent)
                parentRho = self.rep.parent.image(g);
                proj = self.rep.U*self.rep.U';
                rho = self.rep.image(g);
                self.assert(~replab.isNonZeroMatrix(proj*parentRho - parentRho*proj, replab.Settings.doubleEigTol));
                self.R.assertEqv(self.rep.Uinv*parentRho*self.rep.U, rho);
            end
        end
    end
end
