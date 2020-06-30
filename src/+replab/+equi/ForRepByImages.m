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

    methods (Static)

        function e = make(repC, repR, special)
            if ~isa(repR, 'replab.RepByImages') || ~isa(repC, 'replab.RepByImages')
                e = replab.DispatchNext;
                return
            end
            e = replab.equi.ForRepByImages(repC, repR, special);
        end

    end

    methods

        function self = ForRepByImages(repC, repR, special)
            self@replab.Equivariant(repC, repR, special);
        end

        function [X err] = project(self, X)
            X = full(X);
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
            err = replab.equi.errorModel(X);
        end

    end

end
