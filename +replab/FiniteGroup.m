classdef FiniteGroup < replab.FinitelyGeneratedGroup
    
    methods
        
        %ABSTRACT b = contains(self, g)
        % where g can be an element of the group parents
        %
        %ABSTRACT b = knownOrder(self)
        %
        %ABSTRACT o = order(self)
        %
        %ABSTRACT g = sampleUniformly(self)
        %
        %ABSTRACT w = factorization(self, g)
        %
        %ABSTRACT enum = elements(self)

    end
    
    methods
        
        function s = str(self)
            if self.knownOrder
                h = sprintf('%s instance with %d generators of order %s', class(self), self.nGenerators, strtrim(num2str(self.order)));
            else
                h = sprintf('%s instance with %d generators', class(self), self.nGenerators);
            end
            gens = {};
            for i = 1:self.nGenerators
                genName = char('a' + i - 1);
                gens{end+1} = sprintf('  %s: %s', genName, replab.strOf(self.generators{i}));
            end
            s = strjoin({h gens{:}}, newline);
        end
        
    end

    properties (Access = protected)
        randomBag_ = []; % Generator for random elements
    end
    
    methods
        
        function R = randomBag(self)
            if isequal(self.randomBag_, [])
                self.randomBag_ = replab.RandomBag(self.G, self.generators);
            end
            R = self.randomBag_;
        end
        
    end

end
