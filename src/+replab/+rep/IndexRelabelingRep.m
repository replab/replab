classdef IndexRelabelingRep < replab.Rep
% Representation that permutes the indices of a tensor

    properties
        localDimension % dimension of each subsystem in the tensor space
    end

    methods

        function self = IndexRelabelingRep(group, localDimension)
            assert(isa(group, 'replab.PermutationGroup'));
            n = group.domainSize;
            self.group = group;
            self.dimension = localDimension^n;
            self.isUnitary = true;
            self.localDimension = localDimension;
            self.field = 'R';
        end

        function rho = image(self, g)
            n = self.group.domainSize;
            d = self.dimension;
            dims = self.localDimension*ones(1, n);
            I = permute(reshape(1:d, dims), g);
            I = I(:)';
            rho = full(sparse(I, 1:d, ones(1, d), d, d));
        end

    end

end
