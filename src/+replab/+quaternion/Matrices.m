classdef Matrices < replab.Domain
% Describes the vector space of nR x nC quaternion matrices
    
    properties
        nR % integer: Row size
        nC % integer: Column size
    end
    
    methods
        
        function self = Matrices(nR, nC)
            self.nR = nR;
            self.nC = nC;
        end
        
        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d quaternion matrices', self.nR, self.nC);
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            part1 = randn(self.nR, self.nC);
            parti = randn(self.nR, self.nC);
            partj = randn(self.nR, self.nC);
            partk = randn(self.nR, self.nC);
            X = replab.quaternion.H(part1, parti, partj, partk)/2;
        end
        
    end

end
