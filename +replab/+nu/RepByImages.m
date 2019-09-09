classdef RepByImages < replab.nu.Rep
% A finite dimensional not necessarily unitary representation of a finitely generated group
    properties (SetAccess = protected)
        images % Generator images
        inverseImages % Generator inverse images
    end

    methods
        
        function self = RepByImages(group, field, dimension, images, inverseImages)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.images = images;
            self.inverseImages = inverseImages;
        end

        % Rep

        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            word = self.group.factorization(g);
            rho = eye(self.dimension);
            for i = 1:length(word.indices)
                ind = word.indices(i);
                e = word.exponents(i);
                if e > 0
                    ge = self.images{ind}^e;
                else
                    ge = self.inverseImages{ind}^(-e);
                end
                rho = rho * ge;
            end
        end
        
        function rho = inverseImage(self, g)
            word = self.group.factorization(g);
            rho = eye(self.dimension);
            for i = length(word.indices):-1:1
                ind = word.indices(i);
                e = word.exponents(i);
                if e < 0
                    ge = self.images{ind}^(-e);
                else
                    ge = self.inverseImages{ind}^e;
                end
                rho = rho * ge;
            end
        end

    end
    
end
