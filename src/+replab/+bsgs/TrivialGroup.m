classdef TrivialGroup < replab.FiniteGroup
% The trivial group containing the single element `[]`
    methods
        
        function self = TrivialGroup
            self.identity = [];
            self.generators = {};
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
        
    end
    
    methods (Access = protected)
        
        function o = computeOrder(self)
            o = vpi(1);
        end
        
        function e = computeElements(self)
            e = replab.IndexedFamily.lambda(1, @(i) [], @(x) 1);
        end
        
        function d = computeDecomposition(self)
            d = replab.FiniteGroupDecomposition.trivial(group, {[]});
        end
        
    end
    
end
