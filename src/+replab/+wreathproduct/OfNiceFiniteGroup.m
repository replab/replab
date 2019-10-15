classdef OfNiceFiniteGroup < replab.NiceFiniteGroup & replab.semidirectproduct.OfFiniteGroups & replab.wreathproduct.Common
% Wreath product of a permutation group acting on a nice finite group

    methods
        
        function self = OfNiceFiniteGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.semidirectproduct.OfFiniteGroups(phi);
            self.parent = self;
            assert(isa(A, self.requiredType));
            self.n = n;
            self.A = A;
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

        function p = niceMonomorphismImage(self, w)
            p = self.imprimitivePermutation(w, @(a) self.A.niceMonomorphismImage(a));
        end

    end
    
end
