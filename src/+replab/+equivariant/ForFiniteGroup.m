classdef ForFiniteGroup < replab.Equivariant

    methods (Static)

        function e = make(repC, repR, special)
            if ~isa(repC.group, 'replab.FiniteGroup')
                e = replab.DispatchNext;
                return
            end
            e = replab.equivariant.ForFiniteGroup(repC, repR, special);
        end

    end

    methods

        function self = ForFiniteGroup(repC, repR, special)
            self@replab.Equivariant(repC, repR, special);
        end

        function [X err] = project(self, X)
            X = full(X);
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
