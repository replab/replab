classdef PermutationCellAction < replab.Action
    methods
        function self = PermutationCellAction(G, P)
            d = G.domainSize;
            assert(isa(G, 'replab.PermutationGroup'));
            self.G = G;
            self.P = P;
        end
        function C1 = leftAction(self, perm, C)
            C1 = C;
            C1(perm) = C;
        end
        function C1 = rightAction(self, perm, C)
            C1 = C(perm);
        end
    end
end
