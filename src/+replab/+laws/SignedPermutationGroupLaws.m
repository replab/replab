classdef SignedPermutationGroupLaws < replab.laws.FiniteGroupLaws
% Laws obeyed by signed permutations
    properties
        P % (`+replab.Domain`): Signed domain on which those signed permutations act
    end
    methods
        function self = SignedPermutationGroupLaws(T)
            self@replab.laws.FiniteGroupLaws(T);
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
            morphismLaws = replab.laws.GroupMorphismLaws(@(s) replab.SignedPermutation.toPermutation(s), self.T, SymGrp);
        end
        function actionLaws = laws_naturalAction(self)
            actionLaws = self.T.naturalAction.laws;
        end

        function actionLaws = laws_vectorAction(self)
            actionLaws = self.T.vectorAction.laws;
        end

        function actionLaws = laws_matrixAction(self)
            actionLaws = self.T.matrixAction.laws;
        end
        function law_toMatrix_fromMatrix_T(self, t)
            M = replab.SignedPermutation.toMatrix(t);
            t1 = replab.SignedPermutation.fromMatrix(M);
            self.T.assertEqv(t, t1);
        end
    end
end
