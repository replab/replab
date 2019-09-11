classdef RepByImages < replab.Rep
% A finite dimensional unitary representation of a finite group
%
% It works by representing the finite group as a permutation group
% (if it is not already a permutation group), then using a BSGS construction
% that stores the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `replab.bsgs.Chain`).
%
% If the finite group is not a permutation group, a "nice monomorphism"
% in the sense of GAP is used, see:
% https://www.gap-system.org/Manuals/doc/ref/chap40.html#X7FFD731684606BC6)
    properties (SetAccess = protected)
        images % Generator images
        chain % BSGS chain with images
    end
        
    methods
        
        function self = RepByImages(group, field, dimension, images)
        % Constructs a representation from images of group generators
        %
        % Args:
        %   group (instance of `replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   niceMonomorphism: Injective group homomorphism from `self.group` into a permutation group
        %   images (row cell array of orthonormal/unitary matrices): Images of the generators of `group` in the same order
            assert(isa(group, 'replab.NiceFiniteGroup'));
            nG = group.nGenerators;
            assert(length(images) == nG);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
            switch field
              case 'R'
                J = replab.domain.OrthonormalMatrices(dimension);
              case 'C'
                J = replab.domain.UnitaryMatrices(dimension);
              otherwise
                error('Unknown field');
            end
            niceId = group.niceMonomorphism(group.identity);
            n = length(niceId);
            I = zeros(n, nG);
            for i = 1:nG
                I(:,i) = group.niceMonomorphism(group.generator(i));
            end
            self.chain = replab.bsgs.Chain.makeWithImages(n, I, J, images);
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
            img = self.group.niceMonomorphism(g);
            rho = self.chain.image(img);
        end
        
    end

    methods (Static)
        
        function rep1 = fromRep(rep)
            assert(isa(rep.group, 'replab.NiceFiniteGroup'));
            nG = rep.group.nGenerators;
            images = cell(1, nG);
            for i = 1:nG
                images{i} = rep.image(rep.group.generator(i));
            end
            rep1 = replab.RepByImages(rep.group, rep.field, rep.dimension, images);
        end
        
    end
end
