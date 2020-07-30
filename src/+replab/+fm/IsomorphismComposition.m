classdef IsomorphismComposition < replab.Isomorphism & replab.fm.Composition

    methods

        function self = IsomorphismComposition(second, first)
            self@replab.fm.Composition(second, first)
        end

    end

    methods % Implementations

        function s = preimageElement(t)
            s = first.preimageElement(second.preimageElement(t));
        end

    end

    methods (Access = protected) % Implementations

        function I = computeInverse(self)
            I = replab.fm.compose(self.first.inverse, self.second.inverse);
        end

    end

end
