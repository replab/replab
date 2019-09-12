classdef ForFiniteGroup < replab.Equivariant
        
    methods

        function self = ForFiniteGroup(repR, repC)
            replab.Dispatch.assert(isa(repR.group, 'replab.FiniteGroup'));
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
            assert(isa(self.group, 'replab.FiniteGroup'));
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                S = X;
                els = T{i};
                nEls = length(els);
                for j = 2:nEls
                    % TODO: use inverse instead of adjoint 
                    gX = self.repC.action(g, self.repR.action(g, X)')';
                    S = S + gX;
                end
                X = S/nEls;
            end
        end
        
    end
    
end
