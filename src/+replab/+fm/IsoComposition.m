classdef IsoComposition < replab.FiniteIsomorphism

    properties (SetAccess = protected)
        first % (`+replab.FiniteIsomorphism`): First isomorphism
        second % (`+replab.FiniteIsomorphism`): Second isomorphism
    end

    methods

        function self = IsoComposition(second, first)
            self.target = second.target;
            self.source = first.source;
        end

        function t = image(self, s)
            t = self.second.image(self.first.image(s));
        end

        function s = preimage(self, t)
            t = self.first.preimage(self.second.preimage(s));
        end

    end

end
