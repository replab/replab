classdef SignedPermutationsLaws < replab.FiniteGroupLaws
% Laws obeyed by signed permutations

    properties
        P % (`replab.Domain`): Signed domain on which those signed permutations act
    end

    methods

        function self = SignedPermutationsLaws(T)
            self@replab.FiniteGroupLaws(T);
            self.P = replab.domain.signedIntAsDouble(1, T.domainSize);
        end

        function law_toPermutation_fromPermutation_T(self, t)
            p = self.T.toPermutation(t);
            t1 = self.T.fromPermutation(p);
            self.T.assertEqv(t, t1);
        end

        function morphismLaws = laws_toPermutation(self)
            d = self.T.domainSize;
            SymGrp = replab.Permutations(2*d);
            morphismLaws = replab.GroupMorphismLaws(@(s) self.T.toPermutation(s), self.T, SymGrp);
        end

        function actionLaws = laws_naturalAction(self)
            actionLaws = replab.ActionLaws(self.T.naturalAction);
        end

        function actionLaws = laws_vectorAction(self)
            actionLaws = replab.ActionLaws(self.T.vectorAction);
        end

        function actionLaws = laws_matrixAction(self)
            actionLaws = replab.ActionLaws(self.T.matrixAction);
        end

        function law_toMatrix_fromMatrix_T(self, t)
            M = self.T.toMatrix(t);
            t1 = self.T.fromMatrix(M);
            self.T.assertEqv(t, t1);
        end

    end

end
