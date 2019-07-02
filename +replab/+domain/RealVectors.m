classdef RealVectors < replab.Domain
% Describes the vector space of d-dimensional column vectors
    
    properties
        d; % dimension
    end
    
    methods
        
        function self = RealVectors(d)
            self.d = d;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x 1 real vectors', self.d);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            X = randn(self.d, 1);
        end
        
    end

end
