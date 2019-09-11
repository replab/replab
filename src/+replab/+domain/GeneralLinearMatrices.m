classdef GeneralLinearMatrices < replab.Group
% Describes the group of square invertible real or complex matrices
    
    properties
        field % {'R', 'C'}: Underlying field, real (R) or complex (C) matrices
        n % integer: size
    end
    
    properties (Access = protected)
        parent_; % general, not necessarily invertible matrices
    end
    
    methods
        
        function self = GeneralLinearMatrices(field, n)
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
                self.identity = speye(n);
            else
                self.identity = eye(n);
            end
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
