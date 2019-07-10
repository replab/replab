classdef SignedPermutationVectorAction < replab.Action & replab.StrFun
% Describes the action of signed permutations on vectors by
% permuting the coefficients and flipping their signs
    methods
        function self = SignedPermutationVectorAction(G)
            d = G.domainSize;
            assert(isa(G, 'replab.SignedPermutations'));
            desc = sprintf('Action of signed permutations on vectors of %d elements', d);
            self = self@replab.StrFun(desc, desc);
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, 1, 1, G.domainSize);
        end
        function vec = leftAction(self, signedPerm, vec)
            vec(abs(signedPerm)) = vec .* sign(signedPerm(:));
        end
        function vec1 = rightAction(self, vec, signedPerm)
            perm = abs(signedPerm);
            vec1 = vec(perm);
            vec1 = vec1 .* sign(signedPerm(:));
        end
    end
end
