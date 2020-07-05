classdef Intersection

    properties
        base
        lhs
        rhs
    end

    methods

        function self = Intersection(lhsGroup, rhsGroup)
            lhs = lhsGroup.chain.mutableCopy;
            rhs = rhsGroup.chain.mutableCopy;
            rhs.baseChange(lhs.base);
            base = rhs.base;
            if ~isequal(lhs.base, base)
                lhs.baseChange(base);
            end
            self.lhs = lhs;
            self.rhs = rhs;
            self.base = base;
        end

        function [ok, outdata] = test(self, l, g, indata)
            outdata = [];
            prevRhsInv = indata;
            b2 = prevRhsInv(g(self.base(l)));
            i = self.rhs.iDelta(b2, l);
            ok = (i ~= 0);
            if ok
                uinv = self.rhs.Uinv{l}(:,i);
                outdata = uinv(prevRhsInv);
            end
        end

        function s = subgroup(self)
            prop = @(g) self.rhs.contains(g);
            k = length(self.base);
            tests = cell(1, k);
            identity = 1:self.lhs.n;
            for l = 1:k
                tests{l} = @(g, indata) self.test(l, g, indata);
            end
            s = replab.bsgs.subgroupSearch(self.lhs, prop, tests, identity);
        end

    end

end
