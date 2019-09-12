classdef ForRepByImages < replab.Equivariant

    properties (Access = protected)
        decompositionR % transversal images for repR
        decompositionC % transversal images for repC
    end
    
    methods (Access = protected)
       
        function computeImages(self)
            self.decompositionR = self.repR.chain.V;
            if self.repR == self.repC
                self.decompositionC = self.repC.chain.V;
            else
                U = self.repR.chain.U;
                k = length(U);
                decC = cell(1, k);
                for i = 1:k
                    Ri = U{i};
                    nEl = size(Ri, 2);
                    Ci = cell(1, nEl);
                    for j = 1:nEl
                        Ci{j} = self.repC.chain.image(Ri(:,j)');
                    end
                    decC{i} = Ci;
                end
                self.decompositionC = decC;
            end
        end
        
    end
    
    methods

        function self = ForRepByImages(repR, repC)
            replab.Dispatch.assert(isa(repR, 'replab.RepByImages'));
            replab.Dispatch.assert(isa(repC, 'replab.RepByImages'));
            self = self@replab.Equivariant(repR, repC);
        end
        
        function X = project(self, X)
            if isempty(self.decompositionR) || isempty(self.decompositionC)
                self.computeImages;
            end
            R = self.decompositionR;
            C = self.decompositionC;
            for i = length(R):-1:1
                Ri = R{i};
                Ci = C{i};
                nEls = length(Ri);
                S = X;
                for j = 2:nEls
                    gX = Ri{j} * X * Ci{j}';
                    S = S + gX;
                end
                X = S/nEls;
            end
        end
        
    end
    
end
