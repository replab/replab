classdef RepLaws < replab.Laws

    properties
        rep
        G % group of which rep is a representation
          %C % commutant algebra
          %U % Unitary/orthonormal matrices on R or C
          %M % n x n matrices over R, or C
    end

    methods

        function self = RepLaws(rep)
            self.rep = rep;
            d = self.rep.dimension;
            self.G = rep.group;
            % self.C = rep.commutant;
% $$$             self.M = replab.domain.Matrices(rep.field, d, d);
% $$$             if rep.isUnitary
% $$$                 switch rep.field
% $$$                   case 'R'
% $$$                     self.U = replab.OrthogonalGroup(d);
% $$$                   case 'C'
% $$$                     self.U = replab.UnitaryGroup(d);
% $$$                 end
% $$$             else
% $$$                 self.U = replab.GeneralLinearGroup(rep.field, d);
% $$$             end
        end

        function law_identity_(self)
            self.assert(norm(self.rep.image(self.G.identity) - speye(self.rep.dimension), 'fro') <= self.rep.errorBound);
        end

        function law_composition_GG(self, g1, g2)
            if self.rep.isExact
                rho1 = self.rep.image(g1, 'exact');
                rho2 = self.rep.image(g2, 'exact');
                g12 = self.G.compose(g1, g2);
                rho12 = self.rep.image(g12, 'exact');
                rho1rho2 = rho1 * rho2;
                self.assert(all(all(rho12 == rho1rho2)));
            end
            rho1 = self.rep.image(g1, 'double/sparse');
            rho2 = self.rep.image(g2, 'double/sparse');
            g12 = self.G.compose(g1, g2);
            rho12 = self.rep.image(g12, 'double/sparse');
            rho1rho2 = rho1 * rho2;
            c = self.rep.conditionNumberEstimate;
            b = self.rep.errorBound;
            % A1~ * A2~ - A12~ = (A1 + dA1)(A2 + dA2) - A12 - dA12 =~ (A1A2 - A12) + dA1 A2 + A1 dA2 - dA12
            % norm( . , 'fro') <= norm(dA1 A2, 'fro') + norm(A1 dA2, 'fro') + norm(dA12, 'fro')
            tol = 2*c*b + b;
            self.assert(norm(rho1rho2 - rho12, 'fro') <= tol);
        end

        %        function morphismLaws = laws_asGroupHomomorphism(self)
        %            morphismLaws = replab.laws.GroupMorphismLaws(@(g) self.rep.image(g), self.G, self.U);
        %        end

% $$$         function law_commutes_with_commutant_algebra_GC(self, g, c)
% $$$             rho = self.rep.image(g);
% $$$             self.M.assertEqv(rho*c, c*rho);
% $$$         end
% $$$
% $$$         function law_respects_division_algebra_G(self, g)
% $$$             if isequal(self.rep.field, 'R') && isequal(self.rep.isIrreducible, true) && ~isempty(self.rep.frobeniusSchurIndicator)
% $$$                 rho = self.rep.image(g);
% $$$                 switch self.rep.frobeniusSchurIndicator
% $$$                   case 0
% $$$                     rho1 = replab.domain.ComplexTypeMatrices.project(rho);
% $$$                     self.M.assertEqv(rho, rho1);
% $$$                   case -2
% $$$                     rho1 = replab.domain.QuaternionTypeMatrices.project(rho);
% $$$                     self.M.assertEqv(rho, rho1);
% $$$                   case 1
% $$$                     % do nothing
% $$$                   otherwise
% $$$                     error('Wrong Frobenius Schur indicator');
% $$$                 end
% $$$             end
% $$$         end

    end

end
