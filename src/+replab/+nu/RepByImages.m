classdef RepByImages < replab.nu.Rep
% A finite dimensional not necessarily unitary representation of a finitely generated group
    properties (SetAccess = protected)
        images % Generator images
        inverseImages % Inverses of generator images
        chain % BSGS chain with images
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
            for i = 1:group.nGenerators
                imageI = images{i};
                invImageI = inverseImages{i};
                assert(isequal(imageI, [dimension dimension]));
                assert(isequal(invImageI, [dimension dimension]));
                elements{i} = [imageI invImageI];
            end
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
            self.inverseImages = inverseImages;
            J = replab.nu.GeneralLinearMatricesWithInverses(field, dimension);
            niceId = group.niceMonomorphism(group.identity);
            n = length(niceId);
            I = zeros(n, nG);
            for i = 1:nG
                I(:,i) = group.niceMonomorphism(group.generator(i));
            end
            chain = replab.bsgs.Chain(n, J);
            C.insertStrongGenerators(I, elements);
            C.randomizedSchreierSims;
            cut = @(X) X(:, 1:dimension);
            C.mutableMapImages(replab.domain.GeneralLinearMatrices(field, dimension), cut);
            C.makeImmutable;
            self.chain = C;
        end

        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            img = self.group.niceMonomorphism(g);
            rho = self.chain.image(img);
        end
        
        % TODO: speed up inverseImage, by implementing inverseImage in the BSGS chain

    end
    
end
