classdef Laws < replab.Laws
    properties (SetAccess = protected)
        G % replab.semidirectproduct.OfCompactGroups: Semidirect product group
        H % replab.CompactGroup: Group acting
        M % replab.CompactGroup: Group acted upon
    end
    methods
        function self = Laws(G)
            assert(isa(G, 'replab.semidirect.OfCompactGroups'));
            self.G = G;
            self.H = G.H;
            self.M = G.N;
        end
    end
    methods
        function actionLaws = laws_phi(self)
            actionLaws = replab.ActionLaws(self.G.phi);
        end
        function law_compose_elements_HM(self, h, n)
            gh = {h self.M.identity};
            gn = {self.H.identity n};
            g1 = {h n};
            g2 = self.G.compose(gh, gn);
            self.G.assertEqv(g1, g2);
        end
        function law_compose_homomorphism_compatible_HM(self, h, n)
            gh = {h self.M.identity};
            gn = {self.H.identity n};
            % n1 = phi_h(n) = h n h^-1
            gconj = self.G.compose(gh, self.G.composeWithInverse(gn, gh));
            n1 = gconj{2};
            n2 = self.G.phi.leftAction(h, n);
            self.M.assertEqv(n1, n2);
        end
    end

end
