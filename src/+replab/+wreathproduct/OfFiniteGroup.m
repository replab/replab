classdef OfFiniteGroup < replab.semidirectproduct.OfFiniteGroups & replab.wreathproduct.Common
% Wreath product of a permutation group acting on a finite group

    methods
        
        function self = OfFiniteGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.semidirectproduct.OfFiniteGroups(phi);
            assert(isa(A, self.requiredType));
            self.n = n;
            self.A = A;
        end

        function t = requiredType(self)
            t = 'replab.FiniteGroup';
        end
        
    end
    
end
