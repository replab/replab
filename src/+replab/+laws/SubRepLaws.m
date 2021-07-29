classdef SubRepLaws < replab.laws.RepLaws
% Laws for subrepresentations

    methods

        function self = SubRepLaws(rep)
            self = self@replab.laws.RepLaws(rep);
        end

        function law_injection_and_projection_(self)
            if self.rep.isExact
                self.assert(all(all(self.rep.projection('exact') * self.rep.injection('exact') == eye(self.rep.dimension))));
            end
        end

        function law_projector_commutes_G(self, g)
            if self.rep.isExact
                piA = self.rep.projector('exact');
                img = self.rep.parent.image(g, 'exact');
                self.assert(all(all(piA*img == img*piA)));
            end
            piA = self.rep.projector('double/sparse');
            img = self.rep.parent.image(g, 'double/sparse');
            % || rho piA - piA rho ||_F = || rho D - D rho ||_F
            % where D = piA - piE, piA is self.rep.projector and piE is the exact one
            % then || rho D - D rho ||_F <= ||rho||2 ||D||F + ||D||F ||rho||2
            tol = 2*self.rep.projectorErrorBound*self.rep.parent.conditionNumberEstimate;
            self.assertApproxEqual(piA*img, img*piA, tol);
        end

        function law_relation_with_parent_rep_G(self, g)
            if self.rep.isExact
                rep = self.rep.parent.image(g, 'exact');
                sub1 = self.rep.image(g, 'exact');
                sub2 = self.rep.projection('exact') * self.rep.parent.image(g, 'exact') * self.rep.injection('exact');
                self.assert(all(all(sub1 == sub2)));
            end
        end

    end

end
