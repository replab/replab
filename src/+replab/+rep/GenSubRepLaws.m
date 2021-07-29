classdef GenSubRepLaws < replab.Laws

    properties (SetAccess = protected)
        sub % (`+replab.rep.GenSubRep`): Generic subrepresentation
        G % (`+replab.CompactGroup`): Group of which rep is a representation
    end

    methods

        function self = GenSubRepLaws(sub)
            self.sub = sub;
            self.G = sub.group;
        end

        function law_identity_(self)
            rho = self.sub.image(self.G.identity);
            rho1 = speye(self.sub.dimension);
            self.assertApproxEqual(rho, rho1, self.sub.toSubRep.errorBound); % TODO
        end

        function law_composition_GG(self, g1, g2)
            rho1 = self.sub.image(g1);
            rho2 = self.sub.image(g2);
            g12 = self.G.compose(g1, g2);
            rho12 = self.sub.image(g12);
            rho1rho2 = rho1 * rho2;
            c = self.sub.toSubRep.conditionNumberEstimate;
            b = self.sub.toSubRep.errorBound; % TODO
            % A1~ * A2~ - A12~ = (A1 + dA1)(A2 + dA2) - A12 - dA12 =~ (A1A2 - A12) + dA1 A2 + A1 dA2 - dA12
            % norm( . , 'fro') <= norm(dA1 A2, 'fro') + norm(A1 dA2, 'fro') + norm(dA12, 'fro')
            tol = 2*c*b + b;
            self.assertApproxEqual(rho1rho2, rho12, tol);
        end

        function law_unitary_G(self, g)
            if self.sub.isUnitary
                rho = self.sub.image(g);
                rhoI = self.sub.inverseImage(g);
                self.assertApproxEqual(rho, rhoI', 2*self.sub.toSubRep.errorBound); % TODO
            end
        end

        function law_projector_commutes_G(self, g)
            piA = self.sub.projector;
            img = self.sub.parent.image(g, 'double/sparse');
            % || rho piA - piA rho ||_F = || rho D - D rho ||_F
            % where D = piA - piE, piA is self.rep.projector and piE is the exact one
            % then || rho D - D rho ||_F <= ||rho||2 ||D||F + ||D||F ||rho||2
            tol = 2*self.sub.toSubRep.projectorErrorBound*self.sub.parent.conditionNumberEstimate;
            self.assertApproxEqual(piA*img, img*piA, tol*10); % 10 is safety factor
        end

    end

end
