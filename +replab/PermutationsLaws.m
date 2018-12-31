classdef PermutationsLaws < replab.PermutationGroupLaws
    methods
        function self = PermutationsLaws(T)
            self = self@replab.PermutationGroupLaws(T);
        end
        function law_toMatrix_fromMatrix_T(self, t)
            M = self.T.toMatrix(t);
            t1 = self.T.fromMatrix(M);
            self.T.assertEqv(t, t1);
        end
    end
end
