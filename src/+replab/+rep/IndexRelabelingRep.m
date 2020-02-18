classdef IndexRelabelingRep < replab.Rep
% Representation that permutes the indices of a tensor

    properties
        localDimension % dimension of each subsystem in the tensor space
    end

    methods

        function self = IndexRelabelingRep(group, localDimension)
            assert(isa(group, 'replab.PermutationGroup'));
            % own properties
            self.localDimension = localDimension;
            % replab.Rep immutable
            n = group.domainSize;
            self.group = group;
            self.field = 'R';
            self.dimension = localDimension^n;
            % replab.Rep mutable
            self.isUnitary = true;
        end

        function rho = image_internal(self, g)
            n = self.group.domainSize;
            d = self.dimension;
            dims = self.localDimension*ones(1, n);
            I = permute(reshape(1:d, dims), g);
            I = I(:)';
            rho = sparse(I, 1:d, ones(1, d), d, d);
        end

        function rho = imageInverse_internal(self, g)
            rho = self.image_internal(g);
            rho = rho';
        end

    end

end
