classdef DecompositionEquivariant < replab.Equivariant

    properties (Access = protected)
        decompositionR % transversal images for repR
        decompositionC % transversal images for repC
    end
    
    methods (Access = protected)
       
        function computeDecomposition(self)
            if isa(self.repR, 'replab.RepByImages')
                chainR = self.repR.chain;
                self.decompositionR = chainR.V;
                if self.repR == self.repC
                    self.decompositionC = chainC.Vinv;
                else
                    U = chainR.U;
                    k = length(U);
                    decC = cell(1, k);
                    for i = 1:k
                        Ri = U{i};
                        nEl = size(Ri, 2);
                        Ci = cell(1, nEl);
                        if isa(self.repC, 'replab.RepByImages')
                            for j = 1:nEl
                                Ci{j} = self.repC.chain.image(Ri(:,j)');
                            end
                        else
                            for j = 1:nEl
                                Ci{j} = self.repC.image(self.group.niceMonomorphismPreimage(Ri(:,j)');
                            end
                        end
                        decC{i} = Ci;
                    end
                end
            end
        end
        
    end
    
    methods

        function self = DecompositionEquivariant(repR, repC)
            replab.Dispatch.assert(isa(repR.group, 'replab.FiniteGroup'));
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
        % Projects any nR x nC matrix in the equivariant subspace
            assert(isa(self.group, 'replab.FiniteGroup'));
            T = self.group.decomposition.T;
            for i = length(T):-1:1
                X = self.averageOver(X, T{i});
            end
        end
        
    end
    
end
