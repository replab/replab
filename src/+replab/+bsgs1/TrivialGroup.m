classdef TrivialGroup < replab.FiniteGroup
% The trivial group containing the single element `[]`
    methods
        
        function self = TrivialGroup
            self.identity = [];
            self.generators = {};
            self.order = vpi(1);
        end
        
        % Domain
        
        function b = eqv(self, x, y)
            b = true;
        end
        
        function g = sample(self)
            g = [];
        end
        
        % Monoid
        
        function z = compose(self, x, y)
            z = [];
        end
        
        % Group
        
        function xInv = inverse(self, x)
            xInv = [];
        end
        
        % FiniteGroup
        
        function e = elements(self)
            e = replab.Enumerator.lambda(1, @(i) [], @(x) 1);
        end
        
        function d = decomposition(self)
            d = replab.FiniteGroupDecomposition.trivial(group, {[]});
        end
        
    end
    
end
