classdef PermutationGroupLaws < replab.NiceFiniteGroupLaws
    properties (SetAccess = protected)
        P;
    end
    methods
        function self = PermutationGroupLaws(T)
            self@replab.NiceFiniteGroupLaws(T);
            self.P = replab.domain.intAsDouble(1, T.domainSize);
        end
    end
    methods
        function actionLaws = laws_naturalAction(self)
            actionLaws = replab.FaithfulActionLaws(self.T.naturalAction);
        end
        function actionLaws = laws_vectorAction(self)
            actionLaws = replab.ActionLaws(self.T.vectorAction);
        end
        function actionLaws = laws_matrixAction(self)
            actionLaws = replab.ActionLaws(self.T.matrixAction);
        end
        function law_orbits_TP(self, t, p)
            p1 = t(p);
            orbits = self.T.orbits;
            self.assert(self.T.orbits.blockIndex(p) == self.T.orbits.blockIndex(p1));
        end
    end
end
