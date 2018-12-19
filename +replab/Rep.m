classdef Rep < replab.cat.GroupMorphism

    properties
        dimension; % Representation dimension
        group; % Finite group represented
    end
    
    methods
        
        function str(self)
            disp(sprintf('Representation of dimension %d', self.dimension));
        end
        
        function disp(self)
            disp(self.str);
        end

% $$$         function M = sampleGroupElement(self)
% $$$             M = self.image(self.group.randomElement);
% $$$         end
% $$$         
% $$$         function M = sampleGroupAlgebra(self, nSamples, nRounds)
% $$$             if nargin < 3
% $$$                 nRounds = 4;
% $$$             end
% $$$             if nargin < 2
% $$$                 nSamples = 5;
% $$$             end                
% $$$             M = self.targetCat.identity;
% $$$             for i = 1:nRounds
% $$$                 M1 = randn * self.image(self.group.randomElement);
% $$$                 for j = 2:nSamples
% $$$                     M1 = M1 + randn * self.image(self.group.randomElement);
% $$$                 end
% $$$                 M = M * M1;
% $$$             end
% $$$             M = M / sqrt(nSamples^nRounds);
% $$$         end
                
    end
    
end
