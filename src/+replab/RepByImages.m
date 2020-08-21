classdef RepByImages < replab.Rep
% A finite dimensional representation of a finite group
%
% The finite group must have an isomorphism to a permutation group, so that the images of that representation
% can be reprseented using a BSGS construction that storing the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `+replab.+bsgs.ChainWithImages`).

    properties (SetAccess = protected)
        preimages % (cell(1,\*) of `.group` elements): Preimages
        images_internal % (cell(1,\*) of double(\*,\*), may be sparse): Images
        inverseImages_internal % (cell(1,\*) of double(\*,\*), may be sparse): Inverses of images
    end

    methods

        function self = RepByImages(group, field, dimension, preimages, images)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double(\*,\*), may be sparse or symbolic): Images of the preimages
            assert(isa(group, 'replab.FiniteGroup'));
            assert(isa(images, 'cell') && (isempty(images) || isrow(images)));
            assert(length(images) == group.nGenerators);
            knownUnitary = true;
            nG = group.nGenerators;
            images_internal = cell(1, nG);
            inverseImages_internal = cell(1, nG);
            for i = 1:nG
                g = group.generator(i);
                o = group.elementOrder(g);
                img = images{i};
                assert(isequal(size(img), [dimension dimension]));
                inv_img = replab.util.repeatedSquaring(img, o-1, @(x,y) x*y);
                assert(all(all(img * inv_img == eye(dimension))), 'RepByImages can only be used with exact images. Use symbolic arguments.');
                knownUnitary = knownUnitary && all(all(img == ctranspose(inv_img)));
                images_internal{i} = img;
                inverseImages_internal{i} = inv_img;
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
            self.images_internal = images_internal;
            self.inverseImages_internal = inverseImages_internal;
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
            % if any of the images is symbolic, don't use sparse
            useSparse = ~(any(cellfun(@(i) isa(i, 'sym'), [self.images_internal self.inverseImages_internal])));
            if self.isUnitary
                % for unitary/orthogonal representations, we don't need to compute inverses explicitly
                if self.overR
                    target = replab.OrthogonalGroup(d, useSparse);
                else
                    target = replab.UnitaryGroup(d, useSparse);
                end
                symToDouble = replab.Morphism.lambda(target, target, @(X) double(X)); % remove symbolic toolbox stuff
                res = replab.bsgs.ChainWithImages.make(n, target, nicePreimages, self.images_internal, ...
                                                       symToDouble, [], order);
            else
                % for nonunitary/nonorthogonal representations, we store the inverses alongside the matrix, and use a special
                % group structure to perform computations
                target1 = replab.GeneralLinearGroupWithInverses(self.field, self.dimension, useSparse);
                target2 = replab.GeneralLinearGroup(self.field, self.dimension, useSparse);
                cut = replab.Morphism.lambda(target1, target2, @(X) double(X(:, 1:self.dimension)));
                images = cell(1, nG);
                for i = 1:nG
                    images{i} = [self.images_internal{i} self.inverseImages_internal{i}];
                end
                res = replab.bsgs.ChainWithImages.make(n, target1, nicePreimages, images, cut, ...
                                                       [], order);
            end
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'images_internal';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:length(self.images_internal)
                names{1, end+1} = sprintf('images_internal{%d}', i);
                values{1, end+1} = self.images_internal{i};
            end
        end

        % Rep

        function rho = image_internal(self, g)
            perm = self.group.niceMorphism.imageElement(g);
            rho = self.chain.image(perm);
            assert(~isa(rho, 'sym')); % TODO remove
            if isa(rho, 'sym')
                rho = double(rho);
            end
        end

        function rho = inverseImage_internal(self, g)
            perm = self.group.niceMorphism.imageElement(g);
            rho = self.chain.inverseImage(perm);
            assert(~isa(rho, 'sym')); % TODO remove
            if isa(rho, 'sym')
                rho = double(rho);
            end
        end

    end

    methods (Static)

        function rep1 = fromExactRep(rep)
        % Constructs a `.RepByImages` from an existing representation with exact images
            assert(isa(rep.group, 'replab.NiceFiniteGroup'));
            nG = rep.group.nGenerators;
            images = arrayfun(@(i) rep.image_internal(rep.group.generator(i)), 1:nG, 'uniform', 0);
            rep1 = replab.RepByImages(rep.group, rep.field, rep.dimension, rep.group.generators, images);
        end

        function rep = fromImageFunction(group, field, dimension, imageFun)
        % Constructs a RepByImages representation using a given morphism
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group of which to construct a representation
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Dimension of the representation
        %   imageFun (function_handle): Function that returns a matrix for any element of ``G``
        %
        % Returns:
        %   `+replab.RepByImages`: The constructed representation
            assert(isa(group, 'replab.FiniteGroup'), 'The given group must be a FiniteGroup');
            nG = group.nGenerators;
            images = arrayfun(@(i) imageFun(group.generator(i)), 1:nG, 'uniform', 0);
            rep = replab.RepByImages(group, field, dimension, group.generators, images);
        end

    end

end
