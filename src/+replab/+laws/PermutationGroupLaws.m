classdef PermutationGroupLaws < replab.laws.FiniteGroupLaws
    properties (SetAccess = protected)
        P % (`+replab.Domain`): Domain on which those permutations act
    end
    methods
        function self = PermutationGroupLaws(T)
            self@replab.laws.FiniteGroupLaws(T);
            self.P = replab.domain.intAsDouble(1, T.domainSize);
        end
    end
    methods
        function actionLaws = laws_naturalAction(self)
            actionLaws = self.T.naturalAction.laws;
        end
        function actionLaws = laws_vectorAction(self)
            actionLaws = self.T.vectorAction.laws;
        end
        function actionLaws = laws_matrixAction(self)
            actionLaws = self.T.matrixAction.laws;
        end
        function law_orbits_TP(self, t, p)
            p1 = t(p);
            orbits = self.T.orbits;
            self.assert(self.T.orbits.blockIndex(p) == self.T.orbits.blockIndex(p1));
        end
        function law_elements_are_ordered_(self)
            if self.T.order < 200
                E = self.T.elements.toCell;
                M = zeros(length(E), self.T.domainSize);
                for i = 1:length(E)
                    M(i,:) = E{i};
                end
                assertEqual(M, sortrows(M));
            end
        end
        function law_toMatrix_fromMatrix_T(self, t)
            M = replab.Permutation.toMatrix(t);
            t1 = replab.Permutation.fromMatrix(M);
            self.T.assertEqv(t, t1);
        end
    end
end
