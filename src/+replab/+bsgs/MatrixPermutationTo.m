classdef MatrixPermutationTo < replab.bsgs.Backtrack
% Computes the permutation that sends a matrix to another matrix
%
% We find one ``p`` such that ``T == S(inverse(p), inverse(p))`` or ``S == T(p, p)``.

    properties
        S % (double(domainSize,domainSize)): Source matrix
        T % (double(domainSize,domainSize)): Target matrix
    end

    methods

        function self = MatrixPermutationTo(group, S, T, SStabilizer, TStabilizer, debug)
            if nargin < 6 || isempty(debug)
                debug = false;
            end
            if nargin < 5 || isempty(TStabilizer)
                TStabilizer = group.matrixStabilizer(T);
            end
            if nargin < 4 || isempty(SStabilizer)
                SStabilizer = group.matrixStabilizer(S);
            end
            self@replab.bsgs.Backtrack(group, 1:group.domainSize, TStabilizer, SStabilizer, debug);
            self.S = S;
            self.T = T;
        end

        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = false;
            partial = self.base(1:l);
            images = prev(ul(partial));
            if ~isequal(self.S(partial, partial), self.T(images, images))
                return
            end
            if ~isequal(sort(self.S(:,beta)), sort(self.T(:,b)))
                return
            end
            if ~isequal(sort(self.S(beta,:)), sort(self.T(b,:)))
                return
            end
            ok = true;
        end

        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = all(all(self.S == self.T(g,g)));
        end

    end

end
