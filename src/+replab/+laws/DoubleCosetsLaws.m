classdef DoubleCosetsLaws < replab.Laws
% Law checks for double cosets

    properties (SetAccess = protected)
        G % (`+replab.FiniteGroup`): Group
        L % (`+replab.FiniteGroup`): Left subgroup
        R % (`+replab.FiniteGroup`): Right subgroup
        D % (`+replab.DoubleCosets`): Double cosets to test
    end

    methods

        function self = DoubleCosetsLaws(D)
            self.D = D;
            self.L = D.leftSubgroup;
            self.R = D.rightSubgroup;
            self.G = D.group;
        end

        function law_cosetRepresentative_LGR(self, l, g, r)
            rep1 = self.D.cosetRepresentative(g);
            rep2 = self.D.cosetRepresentative(self.G.composeAll({l, g, r}));
            self.G.assertEqv(rep1, rep2);
        end

        % elements calls tranversal

        function laws = laws_elements(self)
            laws = replab.laws.Collection(cellfun(@(dc) dc.laws, self.D.elements, 'uniform', 0));
        end

    end

end
