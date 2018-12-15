classdef PermAsGroup < replab.cat.Group
% Defines the abstract group structure for permutations
% on a given domain size n
    
    properties (SetAccess = private)
        n;
    end
    
    methods
        
        function self = PermAsGroup(n)
            self.n = n;
            self.identity = replab.Perm.identity(self.n);
        end
        
        function xInv = inverse(self, x)
            xInv = replab.Perm.inverse(x);
        end
        
        function z = compose(self, x, y)
            z = replab.Perm.compose(x, y);
        end
        
        function b = isIdentity(self, x)
            b = replab.Perm.isIdentity(x);
        end
        
    end
    
end
