classdef QuaternionVectors < replab.Domain
% Describes the vector space of d x 1 quaternion vectors
    
    properties
        d; % dimension
    end
    
    methods
        
        function self = QuaternionVectors(d)
            self.d = d;
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x 1 quaternion vectors', self.d);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function X = sample(self)
            part1 = randn(self.d, 1);
            parti = randn(self.d, 1);
            partj = randn(self.d, 1);
            partk = randn(self.d, 1);
            X = replab.Quaternion(part1, parti, partj, partk)/2;
        end
        
    end

end
