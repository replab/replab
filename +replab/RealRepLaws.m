classdef RealRepLaws < replab.Laws
    properties
        rep;
        G; % group of which rep is a representation
        C; % centralizer algebra
        R; % GL_d(R)
    end
    methods
        function self = RealRepLaws(rep)
            self.rep = rep;
            d = self.rep.dimension;
            self.G = rep.group;
            self.C = rep.centralizerAlgebra;
            eqvFun = @(X, Y) ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
            self.R = replab.GroupFun('Mat', eqvFun, @() rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
            
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
