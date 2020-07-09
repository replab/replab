classdef RightCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        R % (`.RightCosets`): Right cosets structure
        G % (`.NiceFiniteGroup`): Group
        H % (`.NiceFiniteGroup`): Subgroup
    end
    methods
        function self = RightCosetsLaws(R)
            self.R = R;
            self.G = R.group;
            self.H = R.subgroup;
        end
    end
    methods
        function law_canonical_representative_is_stable_G(self, g)
            r1 = self.R.canonicalRepresentative(g);
            r2 = self.R.canonicalRepresentative(r1);
            assertEqual(r1, r2);
        end
        function law_canonical_representatives_are_unique_HG(self, h, g)
            hg = self.G.compose(h, g);
            assertEqual(self.R.canonicalRepresentative(g), self.R.canonicalRepresentative(hg));
        end
        function law_transversal_size_(self)
            assertEqual(length(self.R.transversal), double(self.G.order/self.H.order));
        end
    end
end
