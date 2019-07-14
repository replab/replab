classdef PermutationNaturalAction < replab.FaithfulAction & replab.StrFun
% Describes the natural action of permutations on their domain
%
% i.e. a permutation described by a vector of images of length n
% acts on integers 1..n
%
% The behavior is undefined for integers outside that range
% (< 1 or > n)
    methods
        function self = PermutationNaturalAction(G)
            assert(isa(G, 'replab.PermutationGroup'));
            desc = sprintf('Natural permutation action on %d elements', G.domainSize);
            self@replab.StrFun(desc, desc);
            self.G = G;
            self.P = replab.domain.intAsDouble(1, G.domainSize);
        end
        function p = findMovedElement(self, perm)
            for i = 1:length(perm)
                if i ~= perm(i)
                    p = i;
                    return
                end
            end
            p = [];            
        end
        function p1 = leftAction(self, g, p)
            p1 = g(p);
        end
        function p1 = rightAction(self, p, g)
        % Finds the inverse image by walking around the cycle
            p1 = p;
            q = g(p);
            while q ~= p
                p1 = q;
                q = g(q);
            end
        end
    end
end
