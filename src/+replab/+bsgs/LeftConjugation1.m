classdef LeftConjugation1 < replab.bsgs.Backtrack1
% Computes a permutation that conjugates a permutation to another permutation
%
% ``t = g s g^-1``
% take ``tc^-1`` in `.tCentralizer` and ``sc`` in `.sCentralizer`
% ``tc^-1 t tc = g sc s sc^-1 g^-1``
% ``t = tc g sc s sc^-1 g^-1 tc^-1``
% thus ``g -> tc g sc``

    properties
        s % (permutation): Source element
        t % (permutation): Target element
    end

    methods

        function self = LeftConjugation1(group, s, t, sCentralizer, tCentralizer, debug)
            if nargin < 6 || isempty(debug)
                debug = false;
            end
            if nargin < 5 || isempty(tCentralizer)
                tCentralizer = group.centralizer(t);
            end
            if nargin < 4 || isempty(sCentralizer)
                sCentralizer = group.centralizer(s);
            end

            sOrbits = replab.bsgs.permutationOrbits(s);
            lengths = cellfun(@length, sOrbits);
            sOrbits = sOrbits(lengths > 1); % filter singletons
            [~, I] = sort(-cellfun(@length, sOrbits));
            sOrbits = sOrbits(I); % largest orbits first
            prescribedBase = [sOrbits{:}];
            self@replab.bsgs.Backtrack1(group, prescribedBase, tCentralizer, sCentralizer, debug);
            self.s = s;
            self.t = t;
        end


        function ok = test(self, l, prevG, ul)
            if l == 1
                ok = true;
            else
                betaPrev = self.partialBase(l-1);
                beta = self.partialBase(l);
                if beta == self.s(betaPrev)
                    % if beta == s(betaPrev)
                    %  g(beta) == g(s(betaPrev))
                    %  g(beta) == t(g(betaPrev))
                    ok = prevG(ul(beta)) == self.t(prevG(betaPrev));
                else
                    ok = true;
                end
            end
        end

        function ok = prop(self, g)
        % Verifies that compose(g, s) == compose(t, g) or t = g xs g^-1
            ok = all(g(self.s) == self.t(g));
        end

    end

end
