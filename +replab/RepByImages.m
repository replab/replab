classdef RepByImages < replab.Rep
% A finite dimensional unitary representation of a finitely generated group
    properties (SetAccess = protected)
        images; % Generator images
    end
        
    methods
        
        function self = RepByImages(group, field, dimension, images)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
        end

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{end+1} = 'images';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:length(self.images)
                names{end+1} = sprintf('images{%d}', i);
                values{end+1} = self.images{i};
            end
        end

        % Rep
        
        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            word = self.group.factorization(g);
            rho = speye(self.dimension); % TODO: support for other fields
            for i = 1:length(word.indices)
                ind = word.indices(i);
                e = word.exponents(i);
                if e > 0
                    ge = self.images{ind}^e;
                else
                    ge = self.images{ind}'^(-e);
                end
                rho = rho * ge;
            end
        end
        
        % Overloading for optimization purposes
        
        %function sub = subRep(self, U)
        %    d = size(U, 1);
        %    newImages = cellfun(@(im) U * im * U', self.images, 'UniformOutput', false);
        %    sub = replab.RepByImages(self.group, self.field, d, newImages);
        %end
        
    end
    
end
