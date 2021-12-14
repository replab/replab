classdef DivisionAlgebraLaws < replab.laws.DomainLaws
% Law checks for division algebras

    methods

        function self = DivisionAlgebraLaws(S)
            self@replab.laws.DomainLaws(S);
        end

        function law_bases_are_dual_(self)
            B = self.S.basis;
            D = self.S.dualBasis;
            for i = 1:self.S.dimension
                for j = 1:self.S.dimension
                    self.assert(trace(D(:,:,i)'*B(:,:,j)) == double(i == j));
                end
            end
        end

        function law_project_S(self, X)
            X1 = self.S.project(X);
            self.assertApproxEqual(X, X1, 0);
        end

        function law_encode_decode_S(self, X)
            X1 = self.S.decode(X);
            X2 = self.S.encode(X1);
        end

    end

end
