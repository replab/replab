classdef GeneralLinearGroup < replab.cat.Group
    
    properties (SetAccess = protected)
        n;
        complex;
        canEqv = false;
        canHash = false;
        canSample = true;
    end
    
    methods
        
        function self = GeneralLinearGroup(n, complex)
            self.n = n;
            self.complex = complex;
            self.identity = eye(n);
        end
        
        function s = str(self)
            s = sprintf('Invertible %d x %d', self.n, self.n);
            if self.complex
                s = [s ' complex double matrices'];
            else
                s = [s ' real double matrices'];
            end
        end
        
        function b = eqv(self, x, y)
            error('Cannot compare floating point matrices');
        end
        
        function h = hash(self, x)
            error('Cannot hash floating point matrices');
        end
        
        function s = sample(self)
            n = self.n;
            if self.complex
                s = randn(n, n) + 1i*randn(n, n);
            else
                s = randn(n, n);
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
