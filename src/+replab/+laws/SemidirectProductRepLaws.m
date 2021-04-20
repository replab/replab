classdef SemidirectProductRepLaws < replab.laws.RepLaws
% Laws for semidirect product representations subrepresentations

    methods

        function self = SemidirectProductRepLaws(rep)
            self = self@replab.laws.RepLaws(rep);
        end

        function law_representations_compatibility_G(self, g)
            h = g{1};
            n = g{2};
            rho1 = self.rep.Nrep.image(self.G.phi.leftAction(h, n));
            rho2 = self.rep.Hrep.image(h) * self.rep.Nrep.image(n) * self.rep.Hrep.inverseImage(h);
            self.assertApproxEqual(rho1, rho2, 10*self.rep.errorBound+1e-15); % TODO: proper error bound
        end

    end

end
