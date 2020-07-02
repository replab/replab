classdef LambdaMorphism < replab.Morphism

    properties (SetAccess = protected)
        imageFun % (function_handle): Function computing images
    end

    methods

        function self = LambdaMorphism(source, target, imageFun)
            self.source = source;
            self.target = target;
            self.imageFun = imageFun;
        end

        function t = image(self, s)
            f = self.imageFun;
            t = f(s);
        end

    end

end
