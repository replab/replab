classdef SemidirectProductGroupLaws < replab.Laws

    properties (SetAccess = protected)
        G % (`+replab.SemidirectProductGroup`): Semidirect product group
        H % (`+replab.CompactGroup`): Group acting
        N % (`+replab.CompactGroup`): Group acted upon
    end

    methods
        function self = SemidirectProductGroupLaws(G)
            assert(isa(G, 'replab.SemidirectProductGroup'));
            self.G = G;
            self.H = G.H;
            self.N = G.N;
        end
    end

    methods

        function actionLaws = laws_phi(self)
            actionLaws = self.G.phi.laws;
        end

        function law_compose_elements_HN(self, h, n)
            gh = {h, self.N.identity};
            gn = {self.H.identity, n};
            g1 = {h, n};
            g2 = self.G.compose(gh, gn);
            self.G.assertEqv(g1, g2);
        end

        function law_compose_homomorphism_compatible_HN(self, h, n)
            gh = {h, self.N.identity};
            gn = {self.H.identity, n};
            % n1 = phi_h(n) = h n h^-1
            gconj = self.G.compose(gh, self.G.composeWithInverse(gn, gh));
            n1 = gconj{2};
            n2 = self.G.phi.leftAction(h, n);
            self.N.assertEqv(n1, n2);
        end

    end

end
