classdef GeneralLinearMatrices < replab.Group
% Describes the group of n x n invertible real or complex matrices
    
    properties
        n % size
        field % {'R', 'C'} real or complex matrices
    end
    
    properties (Access = protected)
        parent_; % complex matrices
    end
    
    methods
        
        function self = GeneralLinearMatrices(n, field)
            self.n = n;
            self.field = field;
            switch field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(n, n);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(n, n);
              otherwise
                error('Unknown field');
            end
            self.identity = eye(n);
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d invertible matrices in %s', self.n, self.n, self.field);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.parent_.sample;
            % a generic gaussian matrix is almost always invertible
        end
        
        % Semigroup/monoid/group
        
        function Z = compose(self, X, Y)
            Z = X * Y;
        end
        
        function XInv = inverse(self, X)
            XInv = inv(X);
        end
        
    end

end
