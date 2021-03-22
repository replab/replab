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
                self.G.assertEqv(T{i}, self.L.cosetRepresentative(T{i}));
            end
        end

        function law_coset_representative_is_stable_G(self, g)
            r1 = self.L.cosetRepresentative(g);
            r2 = self.L.cosetRepresentative(r1);
            self.G.assertEqv(r1, r2);
        end

        function law_factorize_coset_representative_G(self, g)
            c = self.H.leftCoset(g, self.G);
            [l, r] = c.factorizeShortRepresentativeLetters;
            self.G.assertEqv(self.G.imageLetters(l), r);
            assertTrue(c.contains(r));
        end

        function law_coset_representatives_are_unique_GH(self, g, h)
            gh = self.G.compose(g, h);
            self.G.assertEqv(self.L.cosetRepresentative(g), self.L.cosetRepresentative(gh));
        end

        function law_transversal_size_(self)
            assertEqual(length(self.L.transversal), double(self.G.order/self.H.order));
        end

        function morphismLaws = laws_leftAction(self)
            morphismLaws = replab.laws.FiniteMorphismLaws(self.L.leftAction);
        end

    end

end
