classdef CompositionRep < replab.Rep
% Composition of a representation with a morphism

    properties (SetAccess = protected)
        first % (`+replab.Morphism`): Morphism
        second % (`+replab.Rep`): Representation
    end

    methods

        function self = CompositionRep(first, second)
            self.first = first;
            self.second = second;
        end

        % Rep

        function rho = image_internal(self, g)
            rho = self.second.image_internal(self.first.imageElement(g));
        end

        function rho = imageInverse_internal(self, g)
            rho = self.second.imageInverse_internal(self.first.imageElement(g));
        end

    end

end
