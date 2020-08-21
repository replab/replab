classdef LambdaIsomorphism < replab.Isomorphism

    properties (SetAccess = protected)
        preimageElementFun % (function_handle): Element preimage function
        imageElementFun % (function_handle): Element image function
    end

    methods

        function self = Lambda(source, target, preimageElementFun, imageElementFun)
            self.source = source;
            self.target = target;
            self.preimageElementFun = preimageElementFun;
            self.imageElementFun = imageElementFun;
        end

        function s = preimageElement(self, t)
            f = self.preimageElementFun;
            s = f(t);
        end

        function t = imageElement(self, s)
            f = self.imageElementFun;
            t = f(s);
        end

    end

end
