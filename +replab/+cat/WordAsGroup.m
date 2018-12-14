classdef WordAsGroup < replab.cat.Group
% Defines the abstract group structure for associative words in the
% free group
    
    properties (SetAccess = private)
        identity;
        n;
    end
    
    properties (Constant)
        instance = replab.cat.WordAsGroup;
    end
    
    methods
        
        function self = WordAsGroup
            self.identity = replab.Word.identity;
        end
        
        function xInv = inverse(self, x)
            xInv = inv(x);
        end
        
        function z = compose(self, x, y)
            z = x*y;
        end
        
        function b = isIdentity(self, x)
            b = x.length == 0;
        end
        
    end
    
end
