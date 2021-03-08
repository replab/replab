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

    methods % Implementations

        function b = isExact(self)
            b = true;
        end

    end

    methods (Access = protected)

        function e = computeErrorBound(self)
            e = 0;
        end

        function M = matrixRowAction_double_sparse(self, g, M)
            if issparse(M)
                M = self.image_double_sparse(g) * M;
            else
                n = self.group.domainSize;
                d = self.dimension;
                d2 = size(M, 2);
                dims = [self.localDimension*ones(1, n) d2];
                M = reshape(M, dims);
                g = fliplr(n+1-g); % correct for the order of subsystems
                M = reshape(ipermute(reshape(M, dims), [g n+1]), [d d2]);
            end
        end

        function M = matrixColAction_double_sparse(self, g, M)
            if issparse(M)
                gI = self.group.inverse(g);
                M = M * self.image_double_sparse(gI);
            else
                n = self.group.domainSize;
                d = self.dimension;
                d1 = size(M, 1);
                dims = [d1 self.localDimension*ones(1, n)];
                M = reshape(M, dims);
                g = fliplr(n+1-g); % correct for the order of subsystems
                M = reshape(ipermute(reshape(M, dims), [1 g+1]), [d1 d]);
            end
        end

        function rho = image_double_sparse(self, g)
            n = self.group.domainSize;
            d = self.dimension;
            dims = self.localDimension*ones(1, n);
            g = fliplr(n+1-g); % correct for the order of subsystems
            I = permute(reshape(1:d, dims), g);
            I = I(:)';
            rho = sparse(I, 1:d, ones(1, d), d, d);
        end

        function rho = image_exact(self, g)
            rho = replab.cyclotomic.fromDoubles(self.image_double_sparse(g));
        end

    end

end
