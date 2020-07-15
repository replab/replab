classdef IsoComposition < replab.fm.FiniteComposition & replab.FiniteIsomorphism

    methods

        function self = IsoComposition(second, first)
            self@replab.fm.Composition(second, first)
        end

        function s = preimageElement(t)
            s = first.preimageElement(second.preimageElement(t));
        end

        function I = computeInverse(self)
            I = replab.fm.compose(self.first.inverse, self.second.inverse);
        end

    end

end
