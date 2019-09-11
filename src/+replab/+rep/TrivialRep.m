classdef TrivialRep < replab.Rep
% Describes d copies of the real or complex trivial representation of a group

    properties
        identity % Stored copy of the identity matrix
    end
    
    methods
        
        function self = TrivialRep(group, field, dimension)
            assert(isa(group, 'replab.CompactGroup'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            if replab.Settings.useSparse
                self.identity = speye(self.dimension);
            else
                self.identity = eye(self.dimension);
            end
        end
       
        function rho = image(self, g)
            rho = self.identity;
        end
        
        function rho = sample(self)
            rho = self.identity;
        end            
        
        function v = action(self, g, v)
        % do nothing to v
        end
        
        function X = adjointAction(self, g, X)
        % do nothing to X
        end
        
        function complexRep = complexification(self)
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.TrivialRep(self.group, 'C', self.dimension);
        end
        
        function rep = conj(self)
            rep = self;
        end
        
    end
            
end
