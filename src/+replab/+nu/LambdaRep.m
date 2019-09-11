classdef LambdaRep < replab.nu.Rep
% An implementation of a non unitary representation defined by an image function

    properties (SetAccess = protected)
        imageFun % image function
        inverseImageFun % inverse image function
    end
    
    methods
        
        function self = LambdaRep(group, field, dimension, imageFun, inverseImageFun)
        % Constructs a representation from a group's generator images
        %
        % For documentation of the arguments, see `replab.nu.Rep.lambda`.
            assert(isa(group, 'replab.Group'));
            if nargin < 5
                inverseImageFun = [];
            end
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.imageFun = imageFun;
            self.inverseImageFun = inverseImageFun;
        end

        function rho = image(self, g)
            f = self.imageFun;
            rho = f(g);
        end
        
        function rho = inverseImage(self, g)
            f = self.inverseImageFun;
            if isempty(f)
                rho = inverseImage@replab.nu.Rep.inverseImage(g);
            else
                rho = f(g);
            end
        end
        
    end
    
end
