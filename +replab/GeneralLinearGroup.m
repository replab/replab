classdef GeneralLinearGroup < replab.cat.Group
    
    properties (SetAccess = protected)
        n;
        field;
        canEqv = false;
        canHash = false;
        canSample = true;
    end
    
    methods
        
        function self = GeneralLinearGroup(n, field)
            self.n = n;
            self.field = field;
            switch self.field
              case {'R15', 'C15'}
                self.identity = eye(n);
              case {'R7', 'C7'}
                self.identity = eye(n, 'single');
              otherwise
                error(sprintf('Unknown field %s', field));
            end
        end
        
        function s = str(self)
            s = sprintf('Invertible %d x %d matrices in ', self.n, self.n, self.field);
        end
        
        function b = eqv(self, x, y)
            error('Cannot compare floating point matrices');
        end
        
        function h = hash(self, x)
            error('Cannot hash floating point matrices');
        end
        
        function s = sample(self)
            n = self.n;
            switch self.field
              case 'R7'
                s = single(randn(n, n));
              case 'C7'
                s = single(randn(n, n)) + 1i*single(randn(n, n));
              case 'R15'
                s = randn(n, n);
              case 'C15'
                s = randn(n, n) + 1i*randn(n, n);
            end
        end
        
        function z = compose(self, x, y)
            z = x * y;
        end
        
        function y = inverse(self, x)
            y = inv(x);
        end
        
    end
    
end
