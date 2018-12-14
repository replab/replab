classdef MatAsGroup < replab.cat.Group
% Defines the abstract group structure for permutations
% on a given domain size n
    
    properties (SetAccess = private)
        identity;
        n;
    end
    
    methods
        
        function self = MatAsGroup(n)
            self.n = n;
            self.identity = eye(n);
        end
        
        function xInv = inverse(self, x)
            xInv = inv(x);
        end
        
        function z = compose(self, x, y)
            z = x*y;
        end
        
        function b = isIdentity(self, x)
            b = isequal(x, self.identity);
        end
        
    end
    
end
