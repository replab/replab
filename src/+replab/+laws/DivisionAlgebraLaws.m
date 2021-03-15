classdef DivisionAlgebraLaws < replab.laws.DomainLaws
% Law checks for division algebras

    methods

        function self = DivisionAlgebraLaws(T)
            self@replab.laws.DomainLaws(T);
        end

        function law_bases_are_dual_(self)
            B = self.T.basis;
            D = self.T.dualBasis;
            for i = 1:self.T.dimension
                for j = 1:self.T.dimension
                    self.assert(trace(D(:,:,i)'*B(:,:,j)) == double(i == j));
                end
            end
        end

        function law_project_T(self, X)
            X1 = self.T.project(X);
            self.assertApproxEqual(X, X1, 0);
        end

        function law_encode_decode_T(self, X)
            X1 = self.T.decode(X);
            X2 = self.T.encode(X1);
        end

    end

end
