classdef SignedPermutationGroupLaws < replab.FiniteGroupLaws
% Laws obeyed by signed permutations
    properties
        P % (`replab.Domain`): Signed domain on which those signed permutations act
    end
    methods
        function self = SignedPermutationGroupLaws(T)
            self@replab.FiniteGroupLaws(T);
            self.P = replab.domain.signedIntAsDouble(1, T.domainSize);
        end
        function law_toPermutation_fromPermutation_T(self, t)
            p = replab.SignedPermutation.toPermutation(t);
            t1 = replab.SignedPermutation.fromPermutation(p);
            self.T.assertEqv(t, t1);
        end
        function morphismLaws = laws_toPermutation(self)
            d = self.T.domainSize;
            SymGrp = replab.SymmetricGroup(2*d);
            morphismLaws = replab.GroupMorphismLaws(@(s) replab.SignedPermutation.toPermutation(s), self.T, SymGrp);
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
            M = replab.SignedPermutation.toMatrix(t);
            t1 = replab.SignedPermutation.fromMatrix(M);
            self.T.assertEqv(t, t1);
        end
    end
end
