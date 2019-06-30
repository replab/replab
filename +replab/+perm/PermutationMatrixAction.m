classdef PermutationMatrixAction < replab.Action & replab.StrFun
% Describes the action of permutations on square matrices by simultaneous
% permutations of rows and columns
    methods
        function self = PermutationMatrixAction(G)
            d = G.domainSize;
            assert(isa(G, 'replab.PermutationGroup'));
            desc = sprintf('Action of permutations on %d x %d matrices', d, d);
            self = self@replab.StrFun(@(s, mc) desc);
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, d, 1, G.domainSize);
        end
        function M = leftAction(self, perm, M)
            M(perm, perm) = M;
        end
        function M = rightAction(self, M, perm)
            M = M(perm, perm);
        end
    end
end
