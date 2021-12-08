classdef SignedPermutationGroupLaws < replab.laws.FiniteGroupLaws
% Laws obeyed by signed permutations

    properties (SetAccess = protected)
        P % (`+replab.Domain`): Signed domain on which those signed permutations act
    end

    methods

        function self = SignedPermutationGroupLaws(S)
            self@replab.laws.FiniteGroupLaws(S);
            self.P = replab.domain.signedIntAsDouble(1, S.domainSize);
        end

        function law_toPermutation_fromPermutation_S(self, t)
            p = replab.SignedPermutation.toPermutation(t);
            t1 = replab.SignedPermutation.fromPermutation(p);
            self.S.assertEqv(t, t1);
        end

        function morphismLaws = laws_toPermutation(self)
            d = self.S.domainSize;
            SymGrp = replab.S(2*d);
            morphismLaws = replab.laws.GroupMorphismLaws(@(s) replab.SignedPermutation.toPermutation(s), self.S, SymGrp);
        end

        function actionLaws = laws_naturalAction(self)
            actionLaws = self.S.naturalAction.laws;
        end

        function actionLaws = laws_vectorAction(self)
            actionLaws = self.S.vectorAction.laws;
        end

        function actionLaws = laws_matrixAction(self)
            actionLaws = self.S.matrixAction.laws;
        end

        function law_toMatrix_fromMatrix_S(self, t)
            M = replab.SignedPermutation.toMatrix(t);
            t1 = replab.SignedPermutation.fromMatrix(M);
            self.S.assertEqv(t, t1);
        end

    end

end
