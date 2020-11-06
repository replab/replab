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

        function self = RepByImages_monomial(group, field, dimension, preimages, images, imageGroup, imageElements)
        % Constructs a representation from images of group generators
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/sparse double/intval/cyclotomic(\*,\*)): Exact images of the preimages
        %   imageGroup (`+replab.+perm.GeneralizedSymmetricGroup`): Group of generalized permutations
        %   imageElements (cell(1,\*) of `.imageGroup` elements): Generalized permutations corresponding to the images
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.isUnitary = true;
            self.isExact = true;
            self.preimages = preimages;
            self.images_internal = images;
            self.morphism = group.morphismByImages(imageGroup, 'preimages', preimages, 'images', imageElements);
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = image_internal(self, g)
            gp = self.morphism.imageElement(g);
            rho = self.morphism.target.toCyclotomicMatrix(gp);
        end

        function e = computeErrorBound(self)
            if any(self.morphism.target.m == [1 2 4])
                e = 0;
            else
                e = eps(1)*sqrt(2*self.dimension); % max "dimension" non-zero elements can have an error
            end
        end

    end

    methods % Implementations

        function rho = image(self, g, type)
            if nargin < 3
                type = 'double';
            end
            gp = self.morphism.imageElement(g);
            switch type
              case {'cyclotomic', 'native'}
                rho = self.morphism.target.toCyclotomicMatrix(gp);
              case 'intval'
                rho = intval(self.morphism.target.toCyclotomicMatrix(gp));
              case 'double'
                rho = full(self.morphism.target.toSparseMatrix(gp));
              case 'sparse'
                rho = self.morphism.target.toSparseMatrix(gp);
            end
        end

    end

end
