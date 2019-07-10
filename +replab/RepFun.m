classdef RepFun < replab.Rep
% An implementation of a representation defined by an image function
    properties (SetAccess = protected)
        imageFun; % image function
    end
        
    methods
        
        function self = RepFun(group, field, dimension, imageFun)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.Group'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.imageFun = imageFun;
        end
        
    end
    
end
