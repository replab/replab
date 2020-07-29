classdef FiniteIsomorphismComposition < replab.FiniteIsomorphism & replab.fm.FiniteComposition & replab.fm.IsomorphismComposition

    methods

        function self = FiniteIsomorphismComposition(second, first)
            self@replab.fm.FiniteComposition(second, first)
            self@replab.fm.IsomorphismComposition(second, first);
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
