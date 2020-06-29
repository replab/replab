classdef CompositionMorphism < replab.Morphism

    properties (SetAccess = protected)
        first % (`+replab.Morphism`): First morphism
        second % (`+replab.Morphism`): Second morphism
    end

    methods

        function self = CompositionMorphism(first, second)
            self.target = second.target;
            self.source = first.source;
        end

        function t = image(self, s)
            t = self.second.image(self.first.image(s));
        end

    end

end
