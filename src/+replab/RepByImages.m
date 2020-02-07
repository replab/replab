classdef RepByImages < replab.Rep
% A finite dimensional representation of a finite group
%
% It works by representing the finite group as a permutation group (if it is not already a permutation group),
% then using a BSGS construction that stores the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `+replab.+bsgs.Chain`).
%
% If the finite group is not a permutation group, a "nice monomorphism" in the sense of GAP is used, see:
% https://www.gap-system.org/Manuals/doc/ref/chap40.html#X7FFD731684606BC6)

    properties (SetAccess = protected)
        images_internal % Generator images
        inverseImages_internal % Generator inverse images
    end

    properties (Access = protected)
        chain_ % BSGS chain with images
    end

    methods

        %% Own methods

        function self = RepByImages(group, field, dimension, irrepInfo, images, inverseImages)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.NiceFiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   irrepInfo (`+replab.irreducible.Info`): Information about irreducibility
        %   images (cell(1,*) of double(*,*), may be sparse): Images of the generators of ``group`` in the same order
        %   inverseImages (cell(1,*) of double(*,8), may be sparse): Inverse images of the generators
            assert(isa(group, 'replab.NiceFiniteGroup'));
            assert(isa(images, 'cell') && isrow(inverseImages));
            assert(isa(inverseImages, 'cell') && isrow(inverseImages));
            assert(length(inverseImages) == group.nGenerators);
            isUnitary = true;
            for i = 1:self.group.nGenerators
                isUnitary = isUnitary && isequal(images{i}, inverseImages{i}');
                assert(isequal(size(images{i}), [dimension dimension]));
                assert(isequal(size(inverseImages{i}), [dimension dimension]));
            end
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.isUnitary = isUnitary;
            self.irrepInfo = irrepInfo;
            self.images = images;
            self.inverseImages = inverseImages;
        end

        function c = chain(self)
            if isempty(self.chain_)
                if self.isUnitary
                    if self.overR
                        J = replab.OrthogonalGroup(self.dimension);
                    else
                        J = replab.UnitaryGroup(self.dimension);
                    end
                    niceId = self.group.niceMonomorphismImage(self.group.identity);
                    n = length(niceId);
                    nG = self.group.nGenerators;
                    I = zeros(n, nG);
                    for i = 1:nG
                        I(:,i) = self.group.niceMonomorphismImage(self.group.generator(i));
                    end
                    self.chain_ = replab.bsgs.Chain.makeWithImages(n, I, J, self.images);
                else
                    J = replab.GeneralLinearGroupWithInverses(self.field, self.dimension);
                    niceId = self.group.niceMonomorphismImage(self.group.identity);
                    n = length(niceId);
                    nG = self.group.nGenerators;
                    I = zeros(n, nG);
                    elements = cell(1, nG);
                    for i = 1:nG
                        I(:,i) = self.group.niceMonomorphismImage(self.group.generator(i));
                        elements{i} = [self.images{i} self.inverseImages{i}];
                    end
                    C = replab.bsgs.Chain(n, J);
                    C.insertStrongGenerators(I, elements);
                    C.randomizedSchreierSims;
                    cut = @(X) X(:, 1:self.dimension);
                    C.mutableMapImages(replab.GeneralLinearGroup(self.field, self.dimension), cut);
                    C.makeImmutable;
                    self.chain_ = C;
                end
            end
            c = self.chain_;
        end

        %% Str methods

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'images';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:length(self.images)
                names{1, end+1} = sprintf('images{%d}', i);
                values{1, end+1} = self.images{i};
            end
        end

        %% Rep methods

        function rho = image_internal(self, g)
            img = self.group.niceMonomorphismImage(g);
            rho = self.chain.image(img);
        end

        function rho = inverseImage_internal(self, g)
            img = self.group.niceMonomorphismImage(g);
            rho = self.chain.inverseImage(img);
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
                images{i} = rep.image(g);
                inverseImages{i} = rep.inverseImage(g);
            end
            rep1 = replab.RepByImages(rep.group, rep.field, rep.dimension, rep.irrepInfo, images, inverseImages);
        end

    end
end
