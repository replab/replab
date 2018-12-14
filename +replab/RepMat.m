classdef RepMat < handle
    
    properties
        cat;
        group;
        d; % Representation dimension
        field; % either 'R' or 'C'
        images;
    end
    
    methods
        
        function self = RepMat(group, images, d, field)
            self.cat = replab.cat.MatAsGroup(d);
            self.group = group;
            self.images = images;
            self.field = field;
        end
        
        function M = image(self, permutation)
            word = self.group.factorization(permutation);
            M = self.cat.evaluateWord(word, self.images);
        end
        
    end
    
end
