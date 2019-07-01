classdef OrthonormalMatrices < replab.Group
% Describes the group of n x n orthonormal (real) matrices
    
    properties
        n; % size
    end

    properties (Access = protected)
        parent_; % real matrices
    end

    methods
        
        function self = OrthonormalMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.RealMatrices(n, n);
            self.identity = speye(n);
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d orthonormal matrices', self.n, self.n);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            X = self.parent_.sample;
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end
        
        % Semigroup/monoid/group
        
        function Z = compose(self, X, Y)
            Z = X * Y;
        end
        
        function XInv = inverse(self, X)
            XInv = X';
        end
        
    end

end
