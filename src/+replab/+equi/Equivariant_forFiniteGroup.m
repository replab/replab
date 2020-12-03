classdef Equivariant_forFiniteGroup < replab.Equivariant

    methods

        function self = Equivariant_forFiniteGroup(repC, repR, special)
            self@replab.Equivariant(repC, repR, special);
        end

    end

    methods (Access = protected)

        function X = project_exact(self, X)
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                els = T{i};
                nEls = length(els);
                S = X;
                for j = 2:nEls
                    g = els{j};
                    gX = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'exact'), 'exact');
                    S = S + gX;
                end
                X = S/nEls;
            end
        end

        function [X eX] = project_double_sparse(self, X)
            if nargout > 1
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
                    gX = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'double'), 'double');
                    S = S + gX;
                end
                X = S/nEls;
                if nargout > 1
                    eX = nEls*(eR*cC*sX + cR*eC*sX + eX*cR*cC);
                    if eR == 0 && eC == 0
                        eX = eX + norm(eps(X), 'fro');
                    end
                end
            end
        end

    end

end
