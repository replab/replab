classdef RepByImages < replab.nu.Rep
% A finite dimensional not necessarily unitary representation of a finitely generated group
    properties (SetAccess = protected)
        images % Generator images
        inverseImages % Inverses of generator images
    end
    
    properties (Access = protected)
        chain_ % BSGS chain with images
    end

    methods
        
        function self = RepByImages(group, field, dimension, images, inverseImages)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.NiceFiniteGroup'));
            nG = group.nGenerators;
            assert(length(images) == nG);
            assert(length(inverseImages) == nG);
            assert(isa(images, 'cell'));
            assert(isa(inverseImages, 'cell'));
            elements = cell(1, nG);
            for i = 1:nG
                imageI = images{i};
                invImageI = inverseImages{i};
                assert(isequal(size(imageI), [dimension dimension]));
                assert(isequal(size(invImageI), [dimension dimension]));
            end
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
            self.inverseImages = inverseImages;
        end

        function c = chain(self)
            if isempty(self.chain_)
                J = replab.nu.GeneralLinearMatricesWithInverses(self.field, self.dimension);
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
                C.mutableMapImages(replab.domain.GeneralLinearMatrices(self.field, self.dimension), cut);
                C.makeImmutable;
                self.chain_ = C;
            end
            c = self.chain_;
        end
        
        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            img = self.group.niceMonomorphismImage(g);
            rho = self.chain.image(img);
        end
        
        % TODO: speed up inverseImage, by implementing inverseImage in the BSGS chain

    end
    
end
