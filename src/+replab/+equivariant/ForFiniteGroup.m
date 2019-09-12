classdef ForFiniteGroup < replab.Equivariant
        
    methods

        function self = ForFiniteGroup(repR, repC)
            replab.Dispatch.assert(isa(repR.group, 'replab.FiniteGroup'));
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
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
        end
        
    end
    
end
