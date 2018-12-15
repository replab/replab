classdef PermRep < replab.Rep
    
    properties (SetAccess = private)
        permCat;
        genPermImages; % nG x d double array
                       % where nG = self.group.nGenerators
    end
    
    methods
        
        function prettyDisp(self, spaces)
            disp(sprintf('%s Permutation representation of dimension %d', spaces, self.d));
            for i = 1:size(self.genPermImages, 1)
                disp(sprintf('%s %s: %s', spaces, char('a' + i - 1), num2str(self.genPermImages(i,:))));
            end
        end
        
        function self = PermRep(group, genPermImages, field)
            self.group = group;
            self.d = size(genPermImages, 2);
            self.genPermImages = genPermImages;
            self.permCat = replab.cat.PermAsGroup(self.d);
            self.field = field;
        end
        
        function x = permImage(self, permutation)
            word = self.group.factorization(permutation);
            x = self.permCat.identity;
            for i = 1:length(word.indices)
                g = self.genPermImages(word.indices(i), :);
                e = word.exponents(i);
                we = self.permCat.composeN(g, e);
                x = self.permCat.compose(x, we);
            end
        end
        
        function M = image(self, permutation)
            permRep = self.permImage(permutation);
            M = replab.Perm.matrix(permRep);
        end
        
        function R = decompose(self)
            P = replab.prv.Partition.permutationsOrbits(self.genPermImages);
            conjugateBy = zeros(1, self.d);
            conj = [];
            children = cell(1, P.nBlocks);
            shift = 0;
            for i = 1:P.nBlocks
                block = P.block(i);
                bs = length(block);
                where(block) = 1:bs;
                conj = [conj block];
                childPermImages = where(self.genPermImages(:, block));
                children{i} = replab.PermRep(self.group, childPermImages, self.field);
                shift = shift + bs;
            end
            R = replab.PermConjugateRep(conj, replab.DirectSumRep(self.group, children, self.field));
        end
        
    end
        
end
