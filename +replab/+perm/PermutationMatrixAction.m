classdef PermutationMatrixAction < replab.Action
% Describes the action of permutations on square matrices by simultaneous
% permutations of rows and columns
    methods
        function self = PermutationMatrixAction(G)
            assert(isa(G, 'replab.PermutationGroup'));
            d = G.domainSize;
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, d, 1, G.domainSize);
        end
        function str = headerStr(self)
            d = self.G.domainSize;
            str = sprintf('Action of permutations on %d x %d matrices', d, d);
        end
        function M = leftAction(self, perm, M)
            M(perm, perm) = M;
        end
        function M = rightAction(self, M, perm)
            M = M(perm, perm);
        end
    end
end
