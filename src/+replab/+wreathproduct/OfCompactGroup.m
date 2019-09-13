classdef OfCompactGroup < replab.semidirectproduct.OfCompactGroups & replab.wreathproduct.Common
% Wreath product of a permutation group acting on a compact group
    
    methods
        
        function self = OfCompactGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            assert(isa(A, self.requiredType));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.semidirectproduct.OfCompactGroups(phi);
            self.n = n;
            self.A = A;
        end

        function t = requiredType(self)
            t = 'replab.CompactGroup';
        end
        
    end
    
end
