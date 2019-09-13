classdef GeneralLinearMatricesWithInverses < replab.Group
% Describes the group of n x n invertible real or complex matrices
    
    properties
        field % {'R', 'C'} real or complex matrices
        n % size
    end
    
    properties (Access = protected)
        parent_; % complex matrices
    end
    
    methods
        
        function self = GeneralLinearMatricesWithInverses(field, n)
            self.field = field;
            self.n = n;
            switch field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(n, n);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(n, n);
              otherwise
                error('Unknown field');
            end
            if replab.Settings.useSparse
                self.identity = [speye(n) speye(n)];
            else
                self.identity = [eye(n) eye(n)];
            end
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d invertible matrices in %s', self.n, self.n, self.field);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X(:,1:self.n), Y(:,1:self.n));
        end
        
        function X = sample(self)
            X = self.parent_.sample;
            X = [X inv(X)];
            % a generic gaussian matrix is almost always invertible
        end
        
        % Semigroup/monoid/group
        
        function Z = compose(self, X, Y)
            n = self.n;
            Xinv = X(:,n+1:2*n);
            X = X(:,1:n);
            Yinv = Y(:,n+1:2*n);
            Y = Y(:,1:n);
            Z = [X*Y Yinv*Xinv];
        end
        
        function Xinv = inverse(self, X)
            n = self.n;
            Xinv = X(:,n+1:2*n);
            X = X(:,1:n);
            Xinv = [Xinv X];
        end
        
    end

end
