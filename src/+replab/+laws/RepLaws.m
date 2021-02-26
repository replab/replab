classdef RepLaws < replab.Laws

    properties (SetAccess = protected)
        rep % (`+replab.Rep`): Rep being checked
        G % (`+replab.CompactGroup`): Group of which rep is a representation
        C % (`+replab.Equivariant`): Commutant algebra
        M % (`+replab.Domain`): Matrices
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
            rho = self.rep.image(self.G.identity);
            rho1 = speye(self.rep.dimension);
            self.assertApproxEqual(rho, rho1, self.rep.errorBound);
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
            if ~isempty(self.rep.divisionAlgebraName)
                da = replab.DivisionAlgebra(self.rep.divisionAlgebraName);
                rho = self.rep.image(g);
                rho1 = da.project(rho);
                self.assertApproxEqual(rho, rho1, self.rep.errorBound);
            end
        end

        function law_unitary_G(self, g)
            if self.rep.isUnitary
                rho = self.rep.image(g);
                rhoI = self.rep.inverseImage(g);
                self.assertApproxEqual(rho, rhoI', 2*self.rep.errorBound);
            end
        end

        function law_matrixRowAction_GM(self, g, M)
            c = cond(M);
            M1 = self.rep.image(g) * M;
            M2 = self.rep.matrixRowAction(g, M);
            self.assertApproxEqual(M1, M2, c * self.rep.errorBound);
        end

        function law_matrixColAction_GM(self, g, M)
            c = cond(M);
            M1 = M * self.rep.inverseImage(g);
            M2 = self.rep.matrixColAction(g, M);
            self.assertApproxEqual(M1, M2, c * self.rep.errorBound);
        end

        function law_maximalTorusExponents_(self)
            if self.rep.hasMaximalTorusExponents && self.rep.overC % TODO
                [mu, R] = self.rep.group.reconstruction;
                [powers, partition] = self.rep.maximalTorusExponents;
                t = mu.source.sample; % torus element
                rho1 = self.rep.image(mu.imageElement(t));
                rho2 = diag(prod(bsxfun(@power, exp(2i*pi*t), powers), 2));
                tol = 1e-15*sqrt(self.rep.dimension); % assumption: error on each phase is max 1e-15
                self.assertApproxEqual(rho1, rho2, tol);
            end
        end

    end

end
