classdef PermutationGroupLaws < replab.laws.FiniteGroupLaws

    properties (SetAccess = protected)
        P % (`+replab.Domain`): Domain on which those permutations act
    end

    methods

        function self = PermutationGroupLaws(S)
            self@replab.laws.FiniteGroupLaws(S);
            self.P = replab.domain.intAsDouble(1, S.domainSize);
        end

    end

    methods

        function actionLaws = laws_naturalAction(self)
            actionLaws = self.S.naturalAction.laws;
        end

        function actionLaws = laws_vectorAction(self)
            actionLaws = self.S.vectorAction.laws;
        end

        function actionLaws = laws_matrixAction(self)
            actionLaws = self.S.matrixAction.laws;
        end

        function law_orbits_SP(self, t, p)
            p1 = t(p);
            orbits = self.S.orbits;
            self.assert(self.S.orbits.blockIndex(p) == self.S.orbits.blockIndex(p1));
        end

        function law_elements_are_ordered_(self)
            if self.S.order < 200
                E = self.S.elements;
                M = zeros(length(E), self.S.domainSize);
                for i = 1:length(E)
                    M(i,:) = E{i};
                end
                assertEqual(M, sortrows(M));
            end
        end

        function law_toMatrix_fromMatrix_S(self, t)
            M = replab.Permutation.toMatrix(t);
            t1 = replab.Permutation.fromMatrix(M);
            self.S.assertEqv(t, t1);
        end

    end

end
