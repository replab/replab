classdef RepByImages_monomial < replab.RepByImages
% A unitary, monomial, finite dimensional representation of a finite group
%
% Constructs a morphism from the source group to a group of generalized permutations; doesn't
% accelerate computation of images per se, but provides better error bounds, and accelerate
% the computation of equivariant spaces.

    properties (SetAccess = protected)
        morphism % (`+replab.Morphism`): Morphism to a generalized permutation group that provides the images
    end

    methods

        function self = RepByImages_monomial(group, field, dimension, preimages, images, imageGroup, imageElements, varargin)
        % Constructs a representation from images of group generators
        %
        % Additional keyword arguments are passed to the `+replab.Rep` constructor.
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/cyclotomic(\*,\*)): User-provided images of the preimages
        %   imageGroup (`+replab.+perm.GeneralizedSymmetricGroup`): Group of generalized permutations
        %   imageElements (cell(1,\*) of `.imageGroup` elements): Generalized permutations corresponding to the images
            [args, exists, oldValue] = replab.util.keyValuePairsUpdate(varargin, 'isUnitary', true);
            assert(~exists || isempty(oldValue) || isequal(oldValue, true), 'Monomial representations are unitary');
            imagesErrorBound = zeros(1, length(preimages));
            self@replab.RepByImages(group, field, dimension, preimages, images, imagesErrorBound, args{:});
            self.morphism = group.morphismByImages(imageGroup, 'preimages', preimages, 'images', imageElements);
        end

    end

    methods % Implementations

        function b = isExact(self)
            b = true;
        end

        function rho = image(self, g, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            gp = self.morphism.imageElement(g);
            switch type
              case 'exact'
                rho = self.morphism.target.toCyclotomicMatrix(gp);
              case 'double'
                rho = self.morphism.target.toMatrix(gp);
              case 'double/sparse'
                rho = self.morphism.target.toSparseMatrix(gp);
            end
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rep = computeDouble(self)
            exact = replab.rep.RepByImages_exact(self.group, self.field, self.dimension, self.preimages, self.images);
            rep = double(exact);
        end

        function e = computeErrorBound(self)
            if any(self.morphism.target.m == [1 2 4])
                e = 0;
            else
                e = eps(1)*sqrt(2*self.dimension); % max "dimension" non-zero elements can have an error
            end
        end

    end

end
