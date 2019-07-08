classdef QuaternionMatrices < replab.Domain
% Describes the vector space of nR x nC quaternion matrices
    
    properties
        nR; % row size
        nC; % column size
    end
    
    methods
        
        function self = QuaternionMatrices(nR, nC)
            self.nR = nR;
            self.nC = nC;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d quaternion matrices', self.nR, self.nC);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            part1 = randn(self.nR, self.nC);
            parti = randn(self.nR, self.nC);
            partj = randn(self.nR, self.nC);
            partk = randn(self.nR, self.nC);
            X = replab.Quaternion(part1, parti, partj, partk)/2;
        end
        
    end

end
