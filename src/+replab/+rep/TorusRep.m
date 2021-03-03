classdef TorusRep < replab.Rep
% Representation of a `+replab.TorusGroup` using diagonal matrices

    properties (SetAccess = protected)
        torusMap % (integer(\*,\*)): Morphism from the represented torus group to the torus used for the representation
    end

    methods

        function self = TorusRep(group, torusMap)
            d = size(torusMap, 1);
            td = sum(all(torusMap == 0), 2);
            self@replab.Rep(group, 'C', d, 'isUnitary', true, 'trivialDimension', td);
            self.torusMap = torusMap;
        end

    end

    methods (Access = protected) % Implementations

        function e = computeErrorBound(self)
            e = sqrt(self.dimension)*1e-15; % conservative estimate
        end

        function rho = image_double_sparse(self, g)
            rho = sparse(1:self.dimension, 1:self.dimension, exp(2i*pi*mod(self.torusMap*g, 1)));
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false; % Note: redundant, as `.Rep.isExact` is already false
        end

        function p = invariantBlocks(self)
            p = replab.Partition.finest(self.dimension);
        end

        function b = hasMaximalTorusExponents(self)
            b = true;
        end

        function [powers, blockIndex] = maximalTorusExponents(self)
            powers = self.torusMap;
            blockIndex = 1:self.dimension;
        end

    end

end
