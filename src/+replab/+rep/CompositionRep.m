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
            self.group = first.preimageGroup(second.group);
            self.field = second.field;
            self.dimension = second.dimension;
            self.isUnitary = [];
            self.trivialDimension = [];
            self.isIrreducible = [];
            self.frobeniusSchurIndicator = [];
            self.isDivisionAlgebraCanonical = [];
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
