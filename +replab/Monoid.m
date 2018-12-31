classdef Monoid < replab.Semigroup
% Describes a monoid
    properties (SetAccess = protected)
        identity;
    end
    
    methods % Methods with default implementations
        function b = isIdentity(self, x)
        % Returns true if x is the identity, false otherwise
            b = self.eqv(x, self.identity);
        end
        
    end
end
