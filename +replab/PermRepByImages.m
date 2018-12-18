classdef PermRepByImages < replab.Rep
    
    properties
        permImages;
    end
    
    methods
        
        function self = RepByImages(group, images, T)
            assert(length(images) == group.nGenerators);
            self.T = T;
            self.group = group;
            self.images = images;
        end
        
        function rho = image(self, g)
            word = self.group.factorization(g);
            rho = self.T.identity;
            for i = 1:length(word.indices)
                g = self.images{word.indices(i)};
                e = word.exponents(i);
                ge = self.T.composeN(g, e);
                rho = self.T.compose(rho, ge);
            end            
        end
        
        function d = dimension(self)
            d = self.T.n;
        end
        
        function f = field(self)
            f = self.T.field;
        end
        
        function s = str(self)
            s = sprintf('Representation of dimension %d in %s', self.dimension, self.field);
        end
        
        function disp(self)
            disp(self.str);
        end
                
        function M = sampleGroupElement(self)
            M = self.image(self.group.sample);
        end

        function M = sampleGroupAlgebra(self, nSamples, nRounds)
            if nargin < 3
                nRounds = 4;
             end
             if nargin < 2
                 nSamples = 5;
            end                
            M = self.T.identity;
            for i = 1:nRounds
                M1 = randn * self.image(self.group.sample);
                for j = 2:nSamples
                    M1 = M1 + randn * self.image(self.group.sample);
                end
                M = M * M1;
            end
            M = M / sqrt(nSamples^nRounds);
        end
                
    end
    
end
