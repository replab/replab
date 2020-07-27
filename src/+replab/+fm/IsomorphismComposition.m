classdef IsomorphismComposition < replab.fm.Composition & replab.Isomorphism

    methods

        function self = IsomorphismComposition(second, first)
            self@replab.fm.Composition(second, first)
        end

        function s = preimageElement(t)
            s = first.preimageElement(second.preimageElement(t));
        end

    end

    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.fm.compose(self.first.inverse, self.second.inverse);
        end

    end

end
