classdef Intersection < replab.bsgs.Backtrack
% Computes the intersection of groups
%
% See Section 4.6.6 of
% ``D. Holt et al, Handbook of Computational Group Theory (CRC Press, 2005).``
%
% `.base` is a base for both groups
%
% We use the partial base image test on p. 113 of
% ``G. Butler, Fundamental Algorithms for Permutation Groups. vol. 559 (Springer Berlin Heidelberg, 1991).``
% and in Holt, pp. 124-125

    properties
        rhsChain % (.Chain): Stabilizer chain for the second group given
        prevRhsInv % (cell(1,\*) of permutation): Stores the inverses of ``h`` in Holt pp. 124-125
    end

    methods

        function self = Intersection(lhs, rhs, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            base = unique([lhs.lexChain.base rhs.lexChain.base]);
            self@replab.bsgs.Backtrack(lhs, base, knownSubgroup, knownSubgroup, debug);
            self.prevRhsInv = cell(1, length(base)+1);
            self.prevRhsInv{1} = rhs.identity;
            c = rhs.lexChain.mutableCopy;
            c.baseChange(base);
            c.makeImmutable;
            self.rhsChain = c;
        end

        function ok = test(self, l, prev, ul)
            pri = self.prevRhsInv{l};
            b2 = pri(prev(ul(self.base(l))));
            i = self.rhsChain.iDelta(b2, l);
            ok = (i ~= 0);
            if ok
                uinv = self.rhsChain.Uinv{l}(:,i);
                self.prevRhsInv{l+1} = uinv(pri);
            end
        end

        function ok = prop(self, g)
            ok = self.rhsChain.contains(g);
        end

    end

end
