classdef NURepByImages < replab.NURep
% A finite dimensional unitary representation of a finitely generated group
    properties (SetAccess = protected)
        images % Generator images
        invImages % Generator inverse images
    end

    methods
        
        function self = NURepByImages(group, field, dimension, images, invImages)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
            self.invImages = invImages;
        end

        % Rep

        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            word = self.group.factorization(g);
            rho = eye(self.dimension); % TODO: support for other fields
            for i = 1:length(word.indices)
                ind = word.indices(i);
                e = word.exponents(i);
                if e > 0
                    ge = self.images{ind}^e;
                else
                    ge = self.invImages{ind}^(-e);
                end
                rho = rho * ge;
            end
        end

    end
    
end
