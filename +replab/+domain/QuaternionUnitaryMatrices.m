classdef QuaternionUnitaryMatrices < replab.Group
% Describes the group of n x n quaternion unitary matrices
    
    properties
        n; % matrix size
    end
    
    properties (Access = protected)
        parent_;
    end
    
    methods
        
        function self = QuaternionUnitaryMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.QuaternionMatrices(n, n);
            self.identity = eye(n);
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d quaternion unitary matrices', self.n, self.n);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
        % Generates a n x n quaternion unitary matrix
            X = self.parent_.sample';
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end
        
    end

end
