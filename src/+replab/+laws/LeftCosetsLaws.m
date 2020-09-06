classdef LeftCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        L % (`+replab.LeftCosets`): Left cosets structure
        G % (`+replab.FiniteGroup`): Group
        H % (`+replab.FiniteGroup`): Subgroup
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
                assertEqual(T{i}, self.L.cosetRepresentative(T{i}));
            end
        end
        function law_coset_representative_is_stable_G(self, g)
            r1 = self.L.cosetRepresentative(g);
            r2 = self.L.cosetRepresentative(r1);
            assertEqual(r1, r2);
        end

        function law_coset_representatives_are_unique_GH(self, g, h)
            gh = self.G.compose(g, h);
            assertEqual(self.L.cosetRepresentative(g), self.L.cosetRepresentative(gh));
        end
        function law_transversal_size_(self)
            assertEqual(length(self.L.transversal), double(self.G.order/self.H.order));
        end
        function morphismLaws = laws_leftAction(self)
            morphismLaws = replab.laws.FiniteMorphismLaws(self.L.leftAction);
        end
    end

end
