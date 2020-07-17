classdef PermutationTo < replab.bsgs.Backtrack
% Computes the permutation that sends a vector to another vector
%
% We return the set of ``p`` such that ``t == s(inverse(p))`` or ``s == t(p)``.
%

    properties
        s % (double(1,domainSize)): Source vector
        t % (double(1,domainSize)): Target vector
    end

    methods

        function self = PermutationTo(group, s, t, sStabilizer, tStabilizer, debug)
            if nargin < 6 || isempty(debug)
                debug = false;
            end
            if nargin < 5 || isempty(tStabilizer)
                tStabilizer = group.vectorStabilizer(t);
            end
            if nargin < 4 || isempty(sStabilizer)
                sStabilizer = group.vectorStabilizer(s);
            end
            self@replab.bsgs.Backtrack(group, [], tStabilizer, sStabilizer, debug);
            self.s = s;
            self.t = t;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = s(beta) == t(b);
        end


        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = all(self.s == self.t(g));
        end

    end

end
