classdef WreathProductOfCompactGroups < replab.WreathProductGroup & replab.prods.SemidirectProductOfCompactGroups
% Wreath product of a permutation group acting on a compact group

    methods

        function self = WreathProductOfCompactGroups(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self@replab.WreathProductGroup(phi);
            self@replab.prods.SemidirectProductOfCompactGroups(phi);
            self.n = n;
            self.A = A;
        end

    end

end
