classdef RightCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        R % (`+replab.RightCosets`): Right cosets structure
        G % (`+replab.FiniteGroup`): Group
        H % (`+replab.FiniteGroup`): Subgroup
    end

    methods

        function self = RightCosetsLaws(R)
            self.R = R;
            self.G = R.group;
            self.H = R.subgroup;
        end

    end

    methods

        function law_transversal_representatives_are_canonical_(self)
            T = self.R.transversal;
            for i = 1:length(T)
                assertEqual(T{i}, self.R.cosetRepresentative(T{i}));
            end
        end

        function law_canonical_representative_is_stable_G(self, g)
            r1 = self.R.cosetRepresentative(g);
            r2 = self.R.cosetRepresentative(r1);
            assertEqual(r1, r2);
        end

        function law_canonical_representatives_are_unique_HG(self, h, g)
            hg = self.G.compose(h, g);
            assertEqual(self.R.cosetRepresentative(g), self.R.cosetRepresentative(hg));
        end

        function law_transversal_size_(self)
            assertEqual(length(self.R.transversal), double(self.G.order/self.H.order));
        end

    end

end
