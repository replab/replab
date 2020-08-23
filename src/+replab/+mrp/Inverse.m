classdef Inverse < replab.Isomorphism

    properties (SetAccess = protected)
        of % (`+replab.Isomorphism`): Original isomorphism
    end

    methods

        function self = Inverse(of)
            self.of = of;
            self.source = of.target;
            self.target = of.source;
        end

        function f = inverse(self)
            f = self.of;
        end

        function t = imageElement(self, s)
            t = self.of.preimageElement(s);
        end

        function s = preimageElement(self, t)
            s = self.of.imageElement(t);
        end

    end

end
