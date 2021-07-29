classdef IsomorphismComposition < replab.Isomorphism & replab.mrp.Composition

    methods

        function self = IsomorphismComposition(second, first, imageElementFun)
            self@replab.mrp.Composition(second, first, imageElementFun);
        end

    end

    methods % Implementations

        function s = preimageElement(t)
            s = first.preimageElement(second.preimageElement(t));
        end

    end

    methods (Access = protected) % Implementations

        function I = computeInverse(self)
            I = replab.mrp.compose(self.first.inverse, self.second.inverse);
        end

    end

end
