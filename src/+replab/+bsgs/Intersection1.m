classdef Intersection1 < replab.bsgs.Backtrack1
% Computes the intersection of groups
%
% See Section 4.6.6 of
% ``D. Holt et al, Handbook of Computational Group Theory (CRC Press, 2005).``

    properties
        other0
    end

    methods

        function self = Intersection1(group, other, debug)
            if nargin < 3
                debug = false;
            end
            self@replab.bsgs.Backtrack1(group.lexChain, group.lexChain.base, []);
            c = other.lexChain.mutableCopy;
            c.baseChange(group.lexChain.base);
            c.makeImmutable;
            self.other0 = c;
        end

        function ok = test(self, l, prev, ul)
        % unused
        end

        function ok = prop(self, g)
            ok = self.other0.contains(g);
        end

    end

end
