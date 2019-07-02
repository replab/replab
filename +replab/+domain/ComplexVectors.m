classdef ComplexVectors < replab.Domain
% Describes the vector space of d x 1 complex vectors
    
    properties
        d; % dimension
    end
    
    methods
        
        function self = ComplexVectors(d)
            self.d = d;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x 1 complex vectors', self.d, 1);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            realPart = randn(self.d, 1);
            imagPart = randn(self.d, 1);
            X = (realPart + 1i * imagPart)/sqrt(2);
        end
        
    end

end
