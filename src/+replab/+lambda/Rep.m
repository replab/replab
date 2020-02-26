classdef Rep < replab.Rep
% An implementation of a representation defined by image functions

    properties (SetAccess = protected)
        imageFun % function_handle: Image function
        inverseImageFun % function_handle: Inverse image function
    end

    methods

        function self = Rep(group, field, dimension, imageFun, inverseImageFun)
            assert(isa(group, 'replab.Group'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.imageFun = imageFun;
            self.inverseImageFun = inverseImageFun;
        end

        function rho = image_internal(self, g)
            f = self.imageFun;
            rho = f(g);
        end

        function rho = inverseImage_internal(self, g)
            f = self.inverseImageFun;
            if isempty(f)
                gInv = self.group.inverse(g);
                rho = self.image(gInv);
            else
                rho = f(g);
            end
        end

    end

end
