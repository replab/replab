classdef ComplexifiedRep < replab.Rep
% Representation derived by complexfying a real representation

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation being transformed
    end

    methods

        function self = ComplexifiedRep(parent)
            assert(isequal(parent.field, 'R'));
            self.parent = parent;
            self.group = parent.group;
            self.field = 'C';
            self.dimension = parent.dimension;
            self.isUnitary = parent.isUnitary;
            self.irrepInfo = [];
        end

        function s = headerStr(self)
            s = 'Complexification of representation';
        end

        function M = image(self, g)
            M = self.parent.image(g);
        end

        function M = inverseImage(self, g)
            M = self.parent.inverseImage(g);
        end

    end

end
