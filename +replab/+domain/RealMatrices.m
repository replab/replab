classdef RealMatrices < replab.Domain
% Describes the vector space of nR x nC real matrices
    
    properties
        nR; % row size
        nC; % column size
    end
    
    methods
        
        function self = RealMatrices(nR, nC)
            self.nR = nR;
            self.nC = nC;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d real matrices', self.nR, self.nC);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            X = randn(self.nR, self.nC);
        end
        
    end

end
