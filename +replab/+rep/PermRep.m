classdef PermRep < replab.Rep
% A permutation representation of a finite group
    
    properties (SetAccess = protected)
        imageGroup;
        images;
    end

    methods
        
        function self = PermRep(group, field, dimension, images)
            assert(isa(group, 'replab.FiniteGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.imageGroup = replab.Permutations(dimension);
            self.images = images;
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
            word = self.group.factorization(g);
            x = self.imageGroup.identity;
            for i = 1:length(word.indices)
                g = self.images{word.indices(i)};
                e = word.exponents(i);
                ge = self.imageGroup.composeN(g, e);
                x = self.imageGroup.compose(x, ge);
            end
            rho = replab.Permutations.toMatrix(x);
        end
        
    end
    
    methods (Static)
        
        function permRep = fromRep(rep)
            group = rep.group;
            assert(isa(group, 'replab.FiniteGroup'));
            n = group.nGenerators;
            images = cell(1, n);
            for i = 1:n
                g = group.generator(i);
                try
                    perm = replab.Permutations.fromMatrix(rep.image(g));
                    images{i} = perm;
                catch
                    permRep = [];
                    return
                end
            end
            permRep = replab.rep.PermRep(group, rep.field, rep.dimension, images);
        end
        
    end
    
end
