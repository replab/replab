classdef Float64Rep < Rep
    
    properties (SetAccess = private)
        genMatImages; % Cell array of (double) image matrices
    end
    
    methods
        
        function self = MatRep(group, genMatImages, d, field)
            self.group = group;
            self.genMatImages = genMatImages;
            self.d = d;
            self.field = field;
        end
        
        function M = image(self, permutation)
            word = self.group.factorization(permutation);
            M = self.targetCat.evaluateWord(word, self.genMatImages);
        end
        
    end
    
end
