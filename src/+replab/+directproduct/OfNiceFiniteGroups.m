classdef OfNiceFiniteGroups < replab.NiceFiniteGroup & replab.directproduct.OfFiniteGroups

    methods
        
        function self = OfNiceFiniteGroups(factors)
            self = self@replab.directproduct.OfFiniteGroups(factors);
            self.parent = self;
        end
        
        function t = requiredType(self)
            t = 'replab.NiceFiniteGroup';
        end

        %% Domain methods
        
        function b = eqv(self, x, y)
            b = eqv@replab.directproduct.OfGroups(self, x, y);
        end
        
        %% Monoid methods
        
        function z = compose(self, x, y)
            z = compose@replab.directproduct.OfGroups(self, x, y);
        end
 
        %% Group methods
        
        function xInv = inverse(self, x)
            xInv = inverse@replab.directproduct.OfGroups(self, x);
        end

        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = sampleUniformly@replab.directproduct.OfCompactGroups(self); % force method selection
        end

        %% NiceFiniteGroup methods
        
        function p = niceMonomorphismImage(self, g)
            shift = 0;
            p = [];
            for i = 1:self.nFactors
                pf = self.factor(i).niceMonomorphismImage(g{i});
                p = [p pf+shift];
                shift = shift + length(pf);
            end
        end
        
    end

end
