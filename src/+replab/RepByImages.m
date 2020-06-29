classdef RepByImages < replab.Rep
% A finite dimensional representation of a finite group
%
% It works by representing the finite group as a permutation group (if it is not already a permutation group),
% then using a BSGS construction that stores the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `+replab.+bsgs.ChainWithImages`).
%
% If the finite group is not a permutation group, a "nice monomorphism" in the sense of GAP is used, see:
% https://www.gap-system.org/Manuals/doc/ref/chap40.html#X7FFD731684606BC6)

    properties (SetAccess = protected)
        images_internal % (cell(1,\*) of double(\*,\*), may be sparse): Generator images
        inverseImages_internal % (cell(1,\*) of double(\*,\*), may be sparse): Inverses of generator images
    end

    properties (Access = protected)
        chain_ % (`+replab.+bsgs.ChainWithImages`): BSGS chain with images
    end

    methods

        %% Own methods

        function self = RepByImages(group, field, dimension, images, inverseImages)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.NiceFiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   images (cell(1,\*) of double(\*,\*), may be sparse or symbolic): Images of the generators of ``group`` in the same order
        %   inverseImages (cell(1,\*) of double(\*,\*), may be sparse or symbolic): Inverse images of the generators
            assert(isa(group, 'replab.NiceFiniteGroup'));
            assert(isa(images, 'cell') && (isempty(images) || isrow(images)));
            assert(isa(inverseImages, 'cell') && (isempty(inverseImages) || isrow(inverseImages)));
            assert(length(images) == group.nGenerators);
            assert(length(inverseImages) == group.nGenerators);
            knownUnitary = true;
            isInexact = false;
            for i = 1:group.nGenerators
                knownUnitary = knownUnitary && isequal(images{i}, inverseImages{i}');
                assert(isequal(size(images{i}), [dimension dimension]));
                isInexact = isInexact || ~(isa(images{i}, 'sym') || isequal(images{i}, round(images{i})));
                isInexact = isInexact || ~(isa(inverseImages{i}, 'sym') || isequal(inverseImages{i}, round(inverseImages{i})));
                assert(isequal(size(inverseImages{i}), [dimension dimension]));
            end
            if isInexact
                warning('replab:inexactImages', 'Floating point errors accumulate quickly for approximate images. Use symbolic arguments instead.');
            end
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            % replab.Rep mutable
            if knownUnitary
                self.isUnitary = true;
            end
            self.images_internal = images;
            self.inverseImages_internal = inverseImages;
        end

        function c = chain(self)
            if isempty(self.chain_)
                d = self.dimension;
                n = self.group.niceGroup.domainSize;
                order = self.group.order;
                generators = self.group.niceGroup.generators;
                if self.isUnitary
                    if self.overR
                        target = replab.OrthogonalGroup(d);
                    else
                        target = replab.UnitaryGroup(d);
                    end
                    symToDouble = replab.Morphism.lambda(target, target, @(X) double(X)); % remove symbolic toolbox stuff
                    self.chain_ = replab.bsgs.ChainWithImages.make(n, target, generators, self.images_internal, ...
                                                                   symToDouble, [], order);
                else
                    target1 = replab.GeneralLinearGroupWithInverses(self.field, self.dimension);
                    target2 = replab.GeneralLinearGroup(self.field, self.dimension);
                    cut = replab.Morphism.lambda(target1, target2, @(X) double(X(:, 1:self.dimension)));
                    nG = self.group.nGenerators;
                    images = cell(1, nG);
                    for i = 1:nG
                        images{i} = [self.images_internal{i} self.inverseImages_internal{i}];
                    end
                    self.chain_ = replab.bsgs.ChainWithImages.make(n, target1, generators, images, cut, ...
                                                                   [], order);
                end
            end
            c = self.chain_;
        end

        %% Str methods

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

        %% Rep methods

        function rho = image_internal(self, g)
            img = self.group.niceMonomorphismImage(g);
            rho = self.chain.image(img);
            if isa(rho, 'sym')
                rho = double(rho);
            end
        end

        function rho = inverseImage_internal(self, g)
            img = self.group.niceMonomorphismImage(g);
            rho = self.chain.inverseImage(img);
            if isa(rho, 'sym')
                rho = double(rho);
            end
        end

    end

    methods (Static)

        function rep1 = fromRep(rep)
            assert(isa(rep.group, 'replab.NiceFiniteGroup'));
            nG = rep.group.nGenerators;
            images = cell(1, nG);
            inverseImages = cell(1, nG);
            for i = 1:nG
                g = rep.group.generator(i);
                images{i} = rep.image_internal(g);
                inverseImages{i} = rep.inverseImage_internal(g);
            end
            rep1 = replab.RepByImages(rep.group, rep.field, rep.dimension, images, inverseImages);
        end

    end
end
