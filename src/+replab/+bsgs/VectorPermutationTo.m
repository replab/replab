classdef VectorPermutationTo < replab.bsgs.Backtrack
% Computes the permutation from the group that sends a vector to another vector
%
% We find one ``p`` such that ``t == s(inverse(p))`` or ``s == t(p)``.
%
% Example:
% >>> s = [1 3 2 4 5];
% >>> t = [1 4 5 2 3];
% >>> S5 = replab.S(5);
% >>> bktrk = replab.bsgs.VectorPermutationTo(S5, s, t);
% >>> p = bktrk.find
%     1     3     2     4     5
% >>> s = [1 3 3 3 1];
% >>> t = [3 3 3 1 1];
% >>> grp = S5.subgroup({[4 3 2 1 5]});
% >>> bktrk = replab.bsgs.VectorPermutationTo(S5, s, t);
% >>> p = bktrk.find
%     4     3     2     1     5

    properties
        s % (double(1,domainSize)): Source vector
        t % (double(1,domainSize)): Target vector
    end

    methods

        function self = VectorPermutationTo(group, s, t, sStabilizer, tStabilizer, debug)
            if nargin < 6 || isempty(debug)
                debug = false;
            end
            if nargin < 5 || isempty(tStabilizer)
                tStabilizer = group.vectorStabilizer(t);
            end
            if nargin < 4 || isempty(sStabilizer)
                sStabilizer = group.vectorStabilizer(s);
            end
            self@replab.bsgs.Backtrack(group, 1:group.domainSize, tStabilizer, sStabilizer, debug);
            self.s = s;
            self.t = t;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = self.s(beta) == self.t(b);
        end


        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = all(self.s == self.t(g));
        end

    end

end
