classdef IndexRelabelingRep < replab.Rep
% Representation that permutes the indices of a tensor

    properties (SetAccess = protected)
        localDimension % (integer): Dimension of each subsystem in the tensor space
    end

    methods

        function self = IndexRelabelingRep(group, localDimension)
            assert(isa(group, 'replab.PermutationGroup'));
            n = group.domainSize;
            self@replab.Rep(group, 'R', localDimension^n, 'isUnitary', true, 'isIrreducible', localDimension > 1 && n > 0);
            % own properties
            self.localDimension = localDimension;
        end

    end

    methods (Access = protected)

        function b = computeIsUnitary(self)
            b = true;
        end

        function rep = computeDouble(self)
            rep = double(replab.RepByImages.fromExactRep(self));
        end

        function rho = image_double_sparse(self, g)
            n = self.group.domainSize;
            d = self.dimension;
            dims = self.localDimension*ones(1, n);
            I = permute(reshape(1:d, dims), g);
            I = I(:)';
            rho = sparse(I, 1:d, ones(1, d), d, d);
        end

        function rho = image_exact(self, g)
            rho = replab.cyclotomic.fromDoubles(self.image_double(g));
        end

    end

end
