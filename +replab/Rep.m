classdef Rep < handle
    
    properties (SetAccess = protected)
        group;
        d; % Representation dimension
        field; % either 'R' or 'C'
    end
    
    methods
        
        function prettyDisp(self, spaces)
            disp(sprintf('%s Generic representation of dimension %d in %s', spaces, self.d, self.field));
        end
        
        function disp(self)
            self.prettyDisp('');
        end
        
        function C = targetCat(self)
            C = replab.cat.MatAsGroup(self.d);
        end
       
        function M = image(self, permutation)
            error('Not implemented');
        end
        
        function M = sampleGroupElement(self)
            M = self.image(self.group.randomElement);
        end
        
        function M = sampleGroupAlgebra(self, nSamples, nRounds)
            if nargin < 3
                nRounds = 4;
            end
            if nargin < 2
                nSamples = 5;
            end                
            M = self.targetCat.identity;
            for i = 1:nRounds
                M1 = randn * self.image(self.group.randomElement);
                for j = 2:nSamples
                    M1 = M1 + randn * self.image(self.group.randomElement);
                end
                M = M * M1;
            end
            M = M / sqrt(nSamples^nRounds);
        end
                
    end
    
end
