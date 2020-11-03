classdef GeneralizedPermutationNaturalRep < replab.Rep
% Natural representation of (subgroups of) the generalized symmetric group
%
% This serves to unify the handling of permutation representations and signed permutation representations.

    methods

        function self = GeneralizedPermutationNaturalRep(group, field)
            assert(isa(group, 'replab.perm.GeneralizedSymmetricSubgroup'));
            self.group = group;
            self.field = field;
            self.dimension = group.n;
            self.isUnitary = true;
        end

    end

    methods % Implementations

        function rho = image_internal(self, g)
            rho = self.group.toSparseMatrix(g);
        end

    end

end
