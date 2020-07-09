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
        function law_transversal_representatives_are_canonical_(self)
            T = self.L.transversal;
            for i = 1:length(T)
                assertEqual(T{i}, self.L.canonicalRepresentative(T{i}));
            end
        end
        function law_canonical_representative_is_stable_G(self, g)
            r1 = self.L.canonicalRepresentative(g);
            r2 = self.L.canonicalRepresentative(r1);
            assertEqual(r1, r2);
        end

        function law_canonical_representatives_are_unique_GH(self, g, h)
            gh = self.G.compose(g, h);
            assertEqual(self.L.canonicalRepresentative(g), self.L.canonicalRepresentative(gh));
        end
        function law_transversal_size_(self)
            assertEqual(length(self.L.transversal), double(self.G.order/self.H.order));
        end
    end
end
