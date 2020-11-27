classdef Equivariant_forFiniteGroup < replab.Equivariant

    methods

        function self = Equivariant_forFiniteGroup(repC, repR, special)
            self@replab.Equivariant(repC, repR, special);
        end

        function [X eX] = project(self, X, type)
            if nargin < 3
                type = 'double';
            end
            if strcmp(type, 'exact')
                actionType = 'exact';
                if isa(X, 'double')
                    X = replab.cyclotomic.fromDoubles(X);
                end
            else
                actionType = 'double/sparse';
            end
            if nargout > 1 && ~strcmp(type, 'exact')
                sX = replab.numerical.norm2UpperBound(X); % estimate of largest singular value on X
                                                          % we assume this decreases under averaging
                eX = 0; % starting error on X
                eR = self.repR.errorBound; % Frobenius error on repR
                eC = self.repC.errorBound; % Frebenius error on repC
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
            end
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                S = X;
                els = T{i};
                nEls = length(els);
                for j = 2:nEls
                    g = els{j};
                    gX = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, actionType), actionType);
                    S = S + gX;
                end
                X = S/nEls;
                if nargout > 1 && ~strcmp(type, 'exact')
                    eX = nEls*(eR*cC*sX + cR*eC*sX + cR*cC*eX);
                end
            end
            if strcmp(type, 'double')
                X = full(X);
            end
        end

    end

end
