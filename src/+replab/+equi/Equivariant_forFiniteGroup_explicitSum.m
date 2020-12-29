classdef Equivariant_forFiniteGroup_explicitSum < replab.Equivariant
% Equivariant space of representations of a finite group, whose projection is computed using naive summation

    methods

        function self = Equivariant_forFiniteGroup_explicitSum(repR, repC, special)
            self@replab.Equivariant(repR, repC, special);
        end

    end

    methods (Access = protected)

        function X1 = project_exact(self, X)
            E = self.group.elements.toCell;
            g = E{1};
            X1 = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'exact'), 'exact');
            for i = 2:length(E)
                g = E{i};
                X1 = X1 + self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'exact'), 'exact');
            end
            X1 = X1/length(E);
        end

        function [X1, err] = project_double_sparse(self, X)
            if nargout > 1
                sX = replab.numerical.norm2UpperBound(X); % estimate of largest singular value on X
                                                          % we assume this decreases under averaging
                eX = 0; % starting error on X
                eR = self.repR.errorBound; % Frobenius error on repR
                eC = self.repC.errorBound; % Frebenius error on repC
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
            end
            E = self.group.elements.toCell;
            g = E{1};
            X1 = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'double'), 'double');
            for i = 2:length(E)
                g = E{i};
                X1 = X1 + self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X, 'double'), 'double');
            end
            X1 = X1/length(E);
            if nargout > 1
                err = eR*cC*sX + cR*eC*sX + eX*cR*cC; % *nEls/nEls
            end
        end

    end

end
