classdef FiniteIsomorphismComposition < replab.FiniteIsomorphism & replab.mrp.FiniteComposition & replab.mrp.IsomorphismComposition

    methods

        function self = FiniteIsomorphismComposition(second, first, imageElementFun)
            self@replab.mrp.FiniteComposition(second, first, imageElementFun);
            self@replab.mrp.IsomorphismComposition(second, first, imageElementFun);
        end

        function s = preimageElement(self, t)
            s = first.preimageElement(second.preimageElement(t));
        end

        function s = preimageRepresentative(self, t)
            s = preimageRepresentative@replab.mrp.FiniteComposition(self, t);
        end

    end

    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.mrp.compose(self.first.inverse, self.second.inverse);
        end

        function K = computeKernel(self)
            K = computeKernel@replab.mrp.FiniteComposition(self);
        end

        function I = computeImage(self)
            I = self.target;
        end

    end

end
