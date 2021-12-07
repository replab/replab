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

        function l = laws_random_right_coset(self)
            g = self.G.sample;
            rc = self.H.rightCoset(g);
            l = rc.laws;
        end

        function law_right_transversal_representatives_are_canonical_(self)
            T = self.R.transversal;
            for i = 1:length(T)
                self.G.assertEqv(T{i}, self.R.cosetRepresentative(T{i}));
            end
        end

        function law_right_canonical_representative_is_stable_G(self, g)
            r1 = self.R.cosetRepresentative(g);
            r2 = self.R.cosetRepresentative(r1);
            self.G.assertEqv(r1, r2);
        end

        function law_right_canonical_representatives_are_unique_HG(self, h, g)
            hg = self.G.compose(h, g);
            self.G.assertEqv(self.R.cosetRepresentative(g), self.R.cosetRepresentative(hg));
        end

        function law_right_transversal_size_(self)
            assertEqual(length(self.R.transversal), double(self.G.order/self.H.order));
        end

    end

end
