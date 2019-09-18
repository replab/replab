classdef OfNiceFiniteGroups < replab.NiceFiniteGroup & replab.directproduct.OfFiniteGroups
% A direct product of nice finite groups
%    
% In particular, the permutation image of an element of a direct product group
% is simply the concatenation of the permutation images of the factors (which
% are nice finite groups themselves).
    
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
            b = eqv@replab.directproduct.OfCompactGroups(self, x, y);
        end
        
        %% Monoid methods
        
        function z = compose(self, x, y)
            z = compose@replab.directproduct.OfCompactGroups(self, x, y);
        end
 
        %% Group methods
        
        function xInv = inverse(self, x)
            xInv = inverse@replab.directproduct.OfCompactGroups(self, x);
        end

        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = sampleUniformly@replab.directproduct.OfCompactGroups(self); % force method selection
        end

        %% NiceFiniteGroup methods
        
        function p = niceMonomorphismImage(self, g)
            shift = 0;
            p = [];
            % concatenates the permutation images of the factors
            for i = 1:self.nFactors
                pf = self.factor(i).niceMonomorphismImage(g{i});
                p = [p pf+shift];
                shift = shift + length(pf);
            end
        end
        
    end

end
