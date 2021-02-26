classdef TorusMorphismLaws < replab.Laws

    properties (SetAccess = protected)
        morphism % (`+replab.Morphism`): Morphism
        S % (`+replab.CompactGroup`): Source group
        T % (`+replab.CompactGroup`): Target group
        U % (`+replab.TorusGroup`): Source torus group
        V % (`+replab.TorusGroup`): Target torus group
        muS % (`+replab.Morphism`): Morphism from U to S
        muT % (`+replab.Morphism`): Morphism from V to T
    end

    methods

        function self = TorusMorphismLaws(morphism)
            S = morphism.source;
            T = morphism.target;
            muS = S.reconstruction;
            muT = T.reconstruction;
            U = muS.source;
            V = muT.source;
            self.morphism = morphism;
            self.S = S;
            self.T = T;
            self.U = U;
            self.V = V;
            self.muS = muS;
            self.muT = muT;
        end

        function law_tensorMap_U(self, u)
            s = self.muS.imageElement(u);
            t1 = self.morphism.imageElement(s);
            v = self.morphism.imageTorusElement(u);
            t2 = self.muT.imageElement(v);
            self.T.assertEqv(t1, t2);
        end

    end
end
