classdef ForFiniteGroup < replab.Equivariant

    methods

        function self = ForFiniteGroup(repR, repC, special)
            if ~isa(repR.group, 'replab.FiniteGroup')
                error('replab:dispatch:tryNext', 'try next');
            end
            self = self@replab.Equivariant(repR, repC, special);
        end

        function [X err] = project(self, X)
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                S = X;
                els = T{i};
                nEls = length(els);
                for j = 2:nEls
                    g = els{j};
                    gX = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X));
                    S = S + gX;
                end
                X = S/nEls;
            end
            err = replab.equivariant.errorModel(X);
        end

    end

end
