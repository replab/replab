classdef FiniteIsomorphismComposition < replab.FiniteIsomorphism & replab.fm.FiniteComposition & replab.fm.IsomorphismComposition

    methods

        function self = FiniteIsomorphismComposition(second, first)
            self@replab.fm.FiniteComposition(second, first)
            self@replab.fm.IsomorphismComposition(second, first);
        end

        function s = preimageElement(self, t)
            s = first.preimageElement(second.preimageElement(t));
        end

        function s = preimageRepresentative(self, t)
            s = preimageRepresentative@replab.fm.FiniteComposition(self, t);
        end

    end

    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.fm.compose(self.first.inverse, self.second.inverse);
        end

        function K = computeKernel(self)
            K = computeKernel@replab.fm.FiniteComposition(self);
        end

    end

end
