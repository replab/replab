classdef SignedPermutationNaturalAction < replab.Action
% Describes the natural action of signed permutations on their domain
%
% i.e. a permutation described by a vector of images of length n,
% with a possible sign, acts on signed integers {-n,..,-1, 1,..,n}
%
% The behavior is undefined for integers outside that range
    methods
        function self = SignedPermutationNaturalAction(G)
            assert(isa(G, 'replab.signed.PermutationGroup'));
            self.G = G;
            self.P = replab.domain.signedIntAsDouble(1, G.domainSize);
        end
        function str = headerStr(self)
            d = self.G.domainSize;
            str = sprintf('Natural signed permutation action on 2*%d elements', d);
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
