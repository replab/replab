classdef Inverse < replab.FiniteIsomorphism

    properties (SetAccess = protected)
        of % (`+replab.FiniteIsomorphism`): Original isomorphism
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

        function t = image(self, s)
            t = self.of.preimage(s);
        end

        function s = preimage(self, t)
            s = self.of.image(t);
        end

    end

end
