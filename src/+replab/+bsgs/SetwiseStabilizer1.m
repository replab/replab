classdef SetwiseStabilizer1 < replab.bsgs.Backtrack1
% Computes the unordered partition stabilizer of a group

    properties
        set % (integer(1,\*)): Set to stabilize
        mask % (logical(1,domainSize)): Whether the domain point is in the set
    end

    methods

        function self = SetwiseStabilizer1(group, set, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            mask = false(1, group.domainSize);
            mask(set) = true;
            self@replab.bsgs.Backtrack1(group, set, knownSubgroup, knownSubgroup, debug);
            self.set = set;
            self.mask = mask;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = self.mask(beta) == self.mask(b);
        end


        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = all(self.mask == self.mask(g));
        end

    end

end
