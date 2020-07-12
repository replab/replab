classdef Composition < replab.FiniteMorphism

    properties (SetAccess = protected)
        first % (`+replab.FiniteMorphism`): First morphism
        second % (`+replab.FiniteMorphism`): Second morphism
    end

    methods

        function self = CompositionMorphism(second, first)
            self.target = second.target;
            self.source = first.source;
        end

        function t = image(self, s)
            t = self.second.image(self.first.image(s));
        end

    end

end
