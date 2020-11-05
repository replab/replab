classdef ExactRepByImages < replab.RepByImages
% A finite dimensional representation of a finite group
%
% The finite group must have an isomorphism to a permutation group, so that the images of that representation
% can be reprseented using a BSGS construction that storing the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `+replab.+bsgs.ChainWithImages`).

    properties (SetAccess = protected)
        inverseImages % (cell(1,\*) of double/sparse double/cyclotomic(\*,\*)): Image inverses
    end

    methods

        function self = ExactRepByImages(group, field, dimension, preimages, images, inverseImages, knownUnitary, chain)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/sparse double/intval/cyclotomic(\*,\*)): Images of the preimages
        %   inverseImages (cell(1,\*) of double/sparse double/intval/cyclotomic(\*,\*)): Images of the inverses of the preimages
        %   knownUnitary (logical, optional): Whether the representation is known, defaults to false
        %   chain (`+replab.+bsgs.ChainWithImages`, optional): BSGS chain
            if nargin < 7 || isempty(knownUnitary)
                knownUnitary = false;
            end
            if nargin < 8
                chain = [];
            end
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            % replab.Rep mutable
            if knownUnitary
                self.isUnitary = true;
            end
            self.preimages = preimages;
            self.images = images;
            self.inverseImages = inverseImages;
        end

        function c = chain(self)
        % Returns a BSGS chain that computes images of this representation
        %
        % The preimages are permutations.
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: The BSGS chain with matrix images
            c = self.cached('chain', @() self.computeChain);
        end

        function res = computeChain(self)
            m = length(self.preimages);
            d = self.dimension;
            n = self.group.niceMorphism.target.domainSize;
            iso = self.group.niceMorphism;
            nicePreimages = cellfun(@(g) iso.imageElement(g), self.preimages, 'uniform', 0);
            order = self.group.order;
            if self.isUnitary
                % for unitary/orthogonal representations, we don't need to compute inverses explicitly
                if self.overR
                    target = replab.OrthogonalGroup(d, true);
                else
                    target = replab.UnitaryGroup(d, true);
                end
                res = replab.bsgs.ChainWithImages.make(n, target, nicePreimages, self.images, ...
                                                       [], [], order);
            else
                % for nonunitary/nonorthogonal representations, we store the inverses alongside the matrix, and use a special
                % group structure to perform computations
                target1 = replab.GeneralLinearGroupWithInverses(self.field, self.dimension, true);
                target2 = replab.GeneralLinearGroup(self.field, self.dimension, true);
                cut = replab.Morphism.lambda(target1, target2, @(X) double(X(:, 1:self.dimension)));
                images = arrayfun(@(i) [self.images{i} self.inverseImages{i}], 1:m, 'uniform', 0);
                res = replab.bsgs.ChainWithImages.make(n, target1, nicePreimages, images, ...
                                                       [], [], order);
                res = res.mapImages(cut);
            end
        end

    end

    methods % Implementations

        % Rep

        function rho = image_internal(self, g)
            perm = self.group.niceMorphism.imageElement(g);
            rho = self.chain.image(perm);
        end

        function rho = inverseImage_internal(self, g)
            perm = self.group.niceMorphism.imageElement(g);
            rho = self.chain.inverseImage(perm);
        end

    end

end
