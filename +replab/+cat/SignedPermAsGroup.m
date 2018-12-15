classdef SignedPermAsGroup < replab.cat.Group
% Defines the abstract group structure for signed permutations
% on a given domain size n
    
    properties (SetAccess = private)
        n;
    end
    
    methods
        
        function self = SignedPermAsGroup(n)
            self.n = n;
            self.identity = 1:n;
        end
        
        function y = inverse(self, x)
            n = self.n;
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end
        
        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end
        
        function b = isIdentity(self, x)
            b = isequal(x, self.identity);
        end
        
    end
    
end
