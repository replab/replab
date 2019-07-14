classdef SignedPermutationNaturalAction < replab.FaithfulAction & replab.StrFun
% Describes the natural action of signed permutations on their domain
%
% i.e. a permutation described by a vector of images of length n,
% with a possible sign, acts on signed integers {-n,..,-1, 1,..,n}
%
% The behavior is undefined for integers outside that range
    methods
        function self = SignedPermutationNaturalAction(G)
            assert(isa(G, 'replab.SignedPermutations'));
            desc = sprintf('Natural signed permutation action on 2*%d elements', G.domainSize);
            self@replab.StrFun(desc, desc);
            self.G = G;
            self.P = replab.domain.signedIntAsDouble(1, G.domainSize);
        end
        function p = findMovedElement(self, signedPerm)
            for i = 1:length(signedPerm)
                if i ~= signedPerm(i)
                    p = i;
                    return
                end
            end
            p = [];            
        end
        function p1 = leftAction(self, g, p)
            p1 = g(abs(p))*sign(p);
        end
        function p1 = rightAction(self, p, g)
        % Finds the inverse image by walking around the cycle
            p1 = abs(p);
            q = abs(g(abs(p)));
            while q ~= abs(p)
                p1 = q;
                q = abs(g(q));
            end
            p1 = p1*sign(p)*sign(g(p1));
        end
    end
end
