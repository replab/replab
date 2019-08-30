classdef NURep < replab.NURep
% An implementation of a non unitary representation defined by an image function

    properties (SetAccess = protected)
        imageFun % image function
    end
    
    methods
        
        function self = NURep(group, field, dimension, imageFun)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.Group'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.imageFun = imageFun;
        end

        function rho = image(self, g)
            f = self.imageFun;
            rho = f(g);
        end
        
    end
    
end
