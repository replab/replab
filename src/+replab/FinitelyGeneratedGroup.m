classdef FinitelyGeneratedGroup < replab.Group
    
    properties (SetAccess = protected)
        generators; % 1 x nG cell array of generators
    end
    
    methods % Abstract methods
        
        function w = factorization(self, g)
        % Returns the factorization of an element of this group as a word in its generators
            f = self.factorizationFun;
            w = f(g);
        end
        
    end

    methods
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'generators';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nGenerators
                names{1, end+1} = sprintf('generator(%d)', i);
                values{1, end+1} = self.generator(i);
            end
        end
        
        function n = nGenerators(self)
        % Returns the number of generators
            n = length(self.generators);
        end

        function p = generator(self, i) 
        % Returns the i-th generator
            p = self.generators{i};
        end
        
        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th generator of this group
            p = self.inverse(self.generators{i});
        end
        
        function b = isTrivial(self)
        % Returns true if this group is trivial (i.e. has only one element)
            b = self.nGenerators == 0;
        end
        
        function x = evaluateWord(self, word)
        % Evaluates the given word using the generators of this group
            x = self.identity;
            for i = 1:length(word.indices)
                g = self.generators{word.indices(i)};
                e = word.exponents(i);
                we = self.composeN(g, e);
                x = self.compose(x, we);
            end
        end
        
        function f = freeGroup(self)
        % Returns the free group with the same number of generators
        % as this group
            f = replab.FreeGroup(self.nGenerators);
        end
        
        function phi = morphismFromFreeGroup(self)
        % Returns a function handle that evaluates words in the
        % generators of this group
            phi = @(w) self.evaluateWord(w);
        end
        
        function rho = rep(self, field, dimension, images)
        % Constructs a finite dimensional real or complex representation of this group
        %
        %     field: 'R' or 'C' for real or complex
        % dimension: representation dimension
        %    images: 1 x n cell array of matrices providing the images
        %            of the group generators
            rho = replab.RepByImages(self, field, dimension, images);
        end

        function rho = permutationRep(self, dimension, permutations)
        % Returns a real permutation representation of this group
        %
        %    dimension: dimension of the representation
        % permutations: row cell array of images of the generators as permutations of size "dimension"
            S = replab.Permutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, permutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end
        
        function rho = signedPermutationRep(self, dimension, signedPermutations)
        % Returns a real signed permutation representation of this group
        %
        %          dimension: dimension of the representation
        % signedPermutations: row cell array of images of the generators as permutations of size "dimension"
            S = replab.SignedPermutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, signedPermutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end
                
    end
    
end
