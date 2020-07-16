classdef Intersection1 < replab.bsgs.Backtrack1
% Computes the intersection of groups
%
% See Section 4.6.6 of
% ``D. Holt et al, Handbook of Computational Group Theory (CRC Press, 2005).``

    properties
        other
        data
    end

    methods

        function self = Intersection1(group, other, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = group.trivialSubgroup;
            end
            base = unique([group.lexChain.base other.lexChain.base]);
            self@replab.bsgs.Backtrack1(group.lexChain, base, knownSubgroup.lexChain, debug);
            self.data = cell(1, length(base)+1);
            self.data{1} = 1:group.domainSize;
            c = other.lexChain.mutableCopy;
            c.baseChange(base);
            c.makeImmutable;
            self.other = c;
        end

        function ok = test(self, l, prev, ul)
            fprintf('%d ', l);
            prevRhsInv = self.data{l};
            b2 = prevRhsInv(prev(ul(self.base(l))));
            i = self.other.iDelta(b2, l);
            ok = (i ~= 0);
            if ok
                uinv = self.other.Uinv{l}(:,i);
                self.data{l+1} = uinv(prevRhsInv);
            end
        end

        function ok = prop(self, g)
            ok = self.other.contains(g);
        end

    end

end
