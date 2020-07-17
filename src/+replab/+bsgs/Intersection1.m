classdef Intersection1 < replab.bsgs.Backtrack1
% Computes the intersection of groups
%
% See Section 4.6.6 of
% ``D. Holt et al, Handbook of Computational Group Theory (CRC Press, 2005).``
%
% `.base` is a base for both `.group` and `.other`
%
% `.group` is a BSGS chain for the first group given
% `.other` is a BSGS chain for the second group
%
% We use the partial base image test on p. 113 of
% ``G. Butler, Fundamental Algorithms for Permutation Groups. vol. 559 (Springer Berlin Heidelberg, 1991).``
% and in Holt, pp. 124-125

    properties
        rhs
        prevRhsInv % stores the inverses of ``h`` in Holt pp. 124-125
    end

    methods

        function self = Intersection1(lhs, rhs, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            base = unique([lhs.lexChain.base rhs.lexChain.base]);
            self@replab.bsgs.Backtrack1(lhs, base, knownSubgroup, knownSubgroup, debug);
            self.prevRhsInv = cell(1, length(base)+1);
            self.prevRhsInv{1} = rhs.identity;
            c = rhs.lexChain.mutableCopy;
            c.baseChange(base);
            c.makeImmutable;
            self.rhs = c;
        end

        function ok = test(self, l, prev, ul)
            pri = self.prevRhsInv{l};
            b2 = pri(prev(ul(self.base(l))));
            i = self.rhs.iDelta(b2, l);
            ok = (i ~= 0);
            if ok
                uinv = self.rhs.Uinv{l}(:,i);
                self.prevRhsInv{l+1} = uinv(pri);
            end
        end

        function ok = prop(self, g)
            ok = self.rhs.contains(g);
        end

    end

end
