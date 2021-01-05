classdef RepLaws < replab.Laws

    properties
        rep
        G % Group of which rep is a representation
        C % Commutant algebra
        M % Matrices o
    end

    methods

        function self = RepLaws(rep)
            self.rep = rep;
            self.G = rep.group;
            self.C = rep.commutant;
            d = self.rep.dimension;
            self.M = replab.domain.Matrices(rep.field, d, d);
        end

        function law_identity_(self)
            self.assertApproxEqual(self.rep.image(self.G.identity), speye(self.rep.dimension), self.rep.errorBound);
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
            self.assertApproxEqual(rho1rho2, rho12, tol);
        end

        function law_commutes_with_commutant_algebra_GM(self, g, M)
            [X, errX] = self.C.project(M);
            condX = norm(X, 2);
            I = self.rep.image(g);
            errI = self.rep.errorBound;
            condI = self.rep.conditionNumberEstimate;
            % (X + dX)(I + dI) - (I + dI)(X + dX) = XI - IX + dX I + X dI - dX I - dI X + dX dI - dI dX
            % 2 * norm(dX, 'fro') * norm(I, 2) + 2 * norm(dI, 'fro') * norm(X, 2)
            tol = 2*errI*condI + 2*errX*condX;
            self.assertApproxEqual(X * I, I * X, tol);
        end

        function law_respects_division_algebra_G(self, g)
            if self.rep.overR && self.rep.knownIrreducible
                rho = self.rep.image(g);
                switch self.rep.frobeniusSchurIndicator
                    case 0
                      rho1 = replab.domain.ComplexTypeMatrices.project(rho);
                  case -2
                    rho1 = replab.domain.QuaternionTypeMatrices.project(rho, 'group');
                  case 1
                    rho1 = rho;
                    % do nothing
                  otherwise
                    error('Wrong Frobenius Schur indicator');
                end
                self.assertApproxEqual(rho, rho1, self.rep.errorBound);
            end
        end

    end

end
