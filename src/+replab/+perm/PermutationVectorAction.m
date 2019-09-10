classdef PermutationVectorAction < replab.Action
% Describes the action of permutations on vectors by permuting the coefficients
    methods
        function self = PermutationVectorAction(G)
            d = G.domainSize;
            assert(isa(G, 'replab.PermutationGroup'));
            self.G = G;
            self.P = replab.domain.intAsDoubleMatrix(d, 1, 1, G.domainSize);
        end
        function str = headerStr(self)
            d = self.G.domainSize;
            str = sprintf('Action of permutations on vectors of %d elements', d);
        end
        function vec1 = leftAction(self, perm, vec)
            vec1(perm) = vec;
            vec1 = vec1(:);
        end
        function vec1 = rightAction(self, vec, perm)
            vec1 = vec(perm);
        end
    end
end
