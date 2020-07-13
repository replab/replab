classdef Lambda < replab.Morphism

    properties (SetAccess = protected)
        imageElementFun % (function_handle): Element image function
    end

    methods

        function self = Lambda(source, target, imageElementFun)
            self.source = source;
            self.target = target;
            self.imageElementFun = imageElementFun;
        end

        function t = imageElement(self, s)
            f = self.imageElementFun;
            t = f(s);
        end

    end

end
