classdef Equivariant_forFiniteGroup_relativeReynolds < replab.Equivariant
% Equivariant space of representations of a finite group, whose projection is computed using relative Reynolds summation

    methods

        function self = Equivariant_forFiniteGroup_relativeReynolds(repR, repC, special)
            self@replab.Equivariant(repR, repC, special);
        end

    end

    methods % Implementations

        % Equivariant

        function b = isExact(self)
            b = self.repR.isExact && self.repC.isExact;
        end

    end


    methods (Access = protected) % Implementations

        function X = project_exact(self, X)
            T = self.group.setProduct.sets;
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
            switch self.special
              case 'hermitian'
                X = (X + X')/2;
              case 'symmetric'
                X = (X + X.')/2;
            end
        end

        function [X, eX] = project_double_sparse(self, X)
            if nargout > 1
                sX = replab.numerical.norm2UpperBound(X); % estimate of largest singular value on X
                                                          % we assume this decreases under averaging
                eX = 0; % starting error on X
                eR = self.repR.errorBound; % Frobenius error on repR
                eC = self.repC.errorBound; % Frebenius error on repC
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
            end
            T = self.group.setProduct.sets;
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
                end
            end
            switch self.special
              case 'hermitian'
                X = (X + X')/2;
              case 'symmetric'
                X = (X + X.')/2;
            end
        end

    end

end
