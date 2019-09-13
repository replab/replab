classdef OfNiceFiniteGroups < replab.NiceFiniteGroup & replab.semidirectproduct.OfFiniteGroups
% Describes an external semidirect product of nice finite groups
    
    methods
        
        function self = OfNiceFiniteGroups(phi)
            self = self@replab.semidirectproduct.OfFiniteGroups(phi);
            self.parent = self;
        end

        function t = requiredType(self)
            t = 'replab.NiceFiniteGroup';
        end

        %% Domain methods

        function b = eqv(self, x, y)
            b = eqv@replab.semidirectproduct.OfCompactGroups(self, x, y);
        end

        %% Monoid methods

        function z = compose(self, x, y)
            z = compose@replab.semidirectproduct.OfCompactGroups(self, x, y);
        end

        %% Group methods

        function xInv = inverse(self, x)
            xInv = inverse@replab.semidirectproduct.OfCompactGroups(self, x);
        end

        %% CompactGroup methods

        function g = sampleUniformly(self)
            g = sampleUniformly@replab.semidirectproduct.OfCompactGroups(self); % force method selection
        end

        %% NiceFiniteGroup methods

        function p = niceMonomorphismImage(self, g)
            h = g{1};
            n = g{2};
            hp = self.H.niceMonomorphismImage(h);
            np = self.N.niceMonomorphismImage(n);
            p = [hp np+length(hp)];
        end

    end

end
