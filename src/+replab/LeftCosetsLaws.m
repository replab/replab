classdef LeftCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        L % (`.PermutationGroupLeftCosets`): Left cosets structure
        G % (`.PermutationGroup`): Group
        H % (`.PermutationGroup`): Subgroup
    end
    methods
        function self = LeftCosetsLaws(L)
            self.L = L;
            self.G = L.group;
            self.H = L.subgroup;
        end
    end
    methods
        function law_canonical_representatives_are_unique_GH(self, g, h)
            gh = self.G.compose(g, h);
            assertEqual(self.L.canonicalRepresentative(g), self.L.canonicalRepresentative(gh));
        end
        function law_transversal_size_(self)
            assertEqual(length(self.L.transversal), double(self.G.order/self.H.order));
        end
    end
end
