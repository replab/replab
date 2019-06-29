classdef ComplexMatrices < replab.Domain
% Describes the vector space of nR x nC complex matrices
    
    properties
        nR; % row size
        nC; % column size
    end
    
    methods
        
        function self = ComplexMatrices(nR, nC)
            self.nR = nR;
            self.nC = nC;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d complex matrices', self.nR, self.nC);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            realPart = randn(self.nR, self.nC);
            imagPart = randn(self.nR, self.nC);
            X = (realPart + 1i * imagPart)/sqrt(2);
        end
        
    end

end
