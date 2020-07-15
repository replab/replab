classdef Rep < replab.Rep
% An implementation of a representation defined by image functions

    properties (SetAccess = protected)
        image_internalFun % function_handle: Image function
        inverseImage_internalFun % function_handle: Inverse image function
    end

    methods

        function self = Rep(group, field, dimension, image_internalFun, inverseImage_internalFun)
            assert(isa(group, 'replab.Group'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.image_internalFun = image_internalFun;
            self.inverseImage_internalFun = inverseImage_internalFun;
        end

        function rho = image_internal(self, g)
            f = self.image_internalFun;
            rho = f(g);
        end

        function rho = inverseImage_internal(self, g)
            f = self.inverseImage_internalFun;
            rho = f(g);
        end

    end

end
