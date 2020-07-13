classdef Composition < replab.Morphism

    properties (SetAccess = protected)
        first % (`+replab.Morphism`): First morphism
        second % (`+replab.Morphism`): Second morphism
    end

    methods

        function self = Composition(second, first)
            self.target = second.target;
            self.source = first.source;
            self.first = first;
            self.second = second;
        end

        function t = imageElement(self, s)
            t = self.second.imageElement(self.first.imageElement(s));
        end

    end

end
