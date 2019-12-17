classdef SignedPermutationMatrixAction < replab.Action
% Describes the action of signed permutations on square matrices by simultaneous permutations of rows and columns and sign flips
    methods
        function self = SignedPermutationMatrixAction(G)
            assert(isa(G, 'replab.signed.PermutationGroup'));
            d = G.domainSize;
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, d, 1, G.domainSize);
        end
        function str = headerStr(self)
            d = self.G.domainSize;
            str = sprintf('Action of signed permutations on %d x %d matrices', d, d);
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
