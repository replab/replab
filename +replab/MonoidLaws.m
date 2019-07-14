classdef MonoidLaws < replab.SemigroupLaws
    methods
        function self = MonoidLaws(T)
            self@replab.SemigroupLaws(T);
        end
        function law_composeN_zero_T(self, x)
            id1 = self.T.composeN(x, 0);
            id2 = self.T.identity;
            self.T.assertEqv(id1, id2);
        end
        function law_identity(self)
            id = self.T.identity;
            self.assert(self.T.isIdentity(id));
        end
    end
end
