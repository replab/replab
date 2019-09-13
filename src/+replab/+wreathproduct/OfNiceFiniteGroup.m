classdef OfNiceFiniteGroup < replab.semidirectproduct.OfNiceFiniteGroups & replab.wreathproduct.Common
% Wreath product of a permutation group acting on a nice finite group

    methods
        
        function self = OfNiceFiniteGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.semidirectproduct.OfNiceFiniteGroups(phi);
            assert(isa(A, self.requiredType));
            self.n = n;
            self.A = A;
        end

        function t = requiredType(self)
            t = 'replab.NiceFiniteGroup';
        end
        
    end
    
end
