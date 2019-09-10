classdef RepByImages1 < replab.Rep
% A finite dimensional unitary representation of a finitely generated group
    properties (SetAccess = protected)
        niceMonomorphism % Injective group homomorphism from `self.group` to a permutation group
        images % Generator images
        chain % BSGS chain with images
    end
        
    methods
        
        function self = RepByImages1(group, field, dimension, niceMonomorphism, images)
        % Constructs a representation from images of group generators
        %
        % Args:
        %   group (instance of `replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   niceMonomorphism: Injective group homomorphism from `self.group` into a permutation group
        %   images (row cell array of orthonormal/unitary matrices): Images of the generators of `group` in the same order
            assert(isa(group, 'replab.FiniteGroup'));
            nG = group.nGenerators;
            assert(length(images) == nG);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.niceMonomorphism = niceMonomorphism;
            self.images = images;
            switch field
              case 'R'
                J = replab.domain.OrthonormalMatrices(dimension);
              case 'C'
                J = replab.domain.UnitaryMatrices(dimension);
              otherwise
                error('Unknown field');
            end
            niceId = niceMonomorphism(group.identity);
            n = length(niceId);
            I = zeros(n, nG);
            for i = 1:nG
                I(:,i) = niceMonomorphism(group.generator(i));
            end
            self.chain = replab.bsgs1.Chain.makeWithImages(n, I, J, images);
        end

        % Str

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

        % Rep
        
        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            img = self.niceMonomorphism(g);
            rho = self.chain.image(img);
        end
        
    end
    
end
