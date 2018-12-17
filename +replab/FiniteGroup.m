classdef FiniteGroup < replab.FinGenGrp
    
    methods
        
        % Implement
        %
        % function b = contains(self, g)
        % where g can be an element of the group parents
        %
        % function o = order(self)
        %
        % function g = sampleUniformly(self)
        %
        % function w = factorization(self, g)

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
