classdef ComplexifiedRep < replab.Rep
% Representation derived by complexfying a real representation

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation being transformed
    end

    methods

        function self = ComplexifiedRep(parent)
            assert(isequal(parent.field, 'R'));
            % own properties
            self.parent = parent;
            % from replab.Rep, immutable
            self.group = parent.group;
            self.field = 'C';
            self.dimension = parent.dimension;
            % from replab.Rep, mutable
            self.isUnitary = parent.isUnitary;
            self.trivialDimension = parent.trivialDimension;
        end

        function s = headerStr(self)
            s = 'Complexification of representation';
        end

        function M = image_internal(self, g)
            M = self.parent.image_internal(g);
        end

        function M = inverseImage_internal(self, g)
            M = self.parent.inverseImage_internal(g);
        end

    end

end
