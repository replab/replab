classdef SignedPermutationMatrixAction < replab.Action & replab.StrFun
% Describes the action of signed permutations on square matrices by simultaneous
% permutations of rows and columns and sign flips
    methods
        function self = SignedPermutationMatrixAction(G)
            d = G.domainSize;
            assert(isa(G, 'replab.SignedPermutations'));
            desc = sprintf('Action of signed permutations on %d x %d matrices', d, d);
            self = self@replab.StrFun(@(s, mc) desc);
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, d, 1, G.domainSize);
        end
        function M = leftAction(self, signedPerm, M)
            minusSign = find(signedPerm < 0);
            if length(minusSign) > 0
                M(minusSign, :) = -M(minusSign, :);
                M(:, minusSign) = -M(:, minusSign);
            end
            M(abs(signedPerm), abs(signedPerm)) = M;
        end
        function M = rightAction(self, M, signedPerm)
            perm = abs(signedPerm);
            M = M(perm, perm);
            minusSign = find(signedPerm < 0);
            if length(minusSign) > 0
                M(minusSign, :) = -M(minusSign, :);
                M(:, minusSign) = -M(:, minusSign);
            end
        end
    end
end
