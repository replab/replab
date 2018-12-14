classdef FiniteGrp < handle
   
    properties
        cat; % Category describing the elements of this group
    end
    
    methods % ABSTRACT
        
        function n = nGenerators(self) % Returns the number of generators
            error('Not implemented');
        end

        function p = generator(self, i) % Returns the i-th generator
            error('Not implemented');
        end
                
        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th generator of this group
            p = self.cat.inverse(self.generator(i));
        end

        function gcell = generators(self)
        % Returns a row cell vector of generators
            gcell = cell(1, self.nGenerators);
            for i = 1:self.nGenerators
                gcell{i} = self.generator(i);
            end
        end
        
        function b = isTrivial(self)
            b = self.nGenerators == 0;
        end

        function p = randomElement(self) 
        % Returns a random element from the group
            p = self.randomBag.sample;
        end

    end
    
    properties (Access = protected)
        randomBag_ = []; % Generator for random elements
    end
    
    methods (Access = protected)
        
        function R = randomBag(self)
            if isequal(self.randomBag_, [])
                self.randomBag_ = replab.prv.RandomBag(self.generators, self.cat);
            end
            R = self.randomBag_;
        end
        
    end
    
end
