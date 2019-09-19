classdef ForRepByImages < replab.Equivariant

    properties (Access = protected)
        decompositionR % transversal images for repR
        decompositionC % transversal inverse images for repC
    end
    
    methods (Access = protected)
       
        function computeImages(self)
            self.decompositionR = self.repR.chain.V;
            if self.repR == self.repC
                self.decompositionC = self.repC.chain.Vinv;
            else
                U = self.repR.chain.U;
                k = length(U);
                decC = cell(1, k);
                for i = 1:k
                    Ri = U{i};
                    nEl = size(Ri, 2);
                    Ci = cell(1, nEl);
                    for j = 1:nEl
                        Ci{j} = self.repC.chain.inverseImage(Ri(:,j)');
                    end
                    decC{i} = Ci;
                end
                self.decompositionC = decC;
            end
        end
        
    end
    
    methods

        function self = ForRepByImages(repR, repC)
            if ~isa(repR, 'replab.RepByImages')
                error('replab:dispatch:tryNext', 'try next');
            end
            if ~isa(repC, 'replab.RepByImages')
                error('replab:dispatch:tryNext', 'try next');
            end
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
            if isempty(self.decompositionR) || isempty(self.decompositionC)
                self.computeImages;
            end
            R = self.decompositionR;
            Cinv = self.decompositionC;
            for i = length(R):-1:1
                Ri = R{i};
                Cinvi = Cinv{i};
                nEls = length(Ri);
                S = X;
                for j = 2:nEls
                    gX = Ri{j} * X * Cinvi{j};
                    S = S + gX;
                end
                X = S/nEls;
            end
        end
        
    end
    
end
