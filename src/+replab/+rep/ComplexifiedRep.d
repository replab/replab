classdef ComplexifiedRep < replab.Rep
% Representation derived by complexfying a real representation

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Real representation being transformed
    end

    methods

        function self = ComplexifiedRep(parent)
            assert(isequal(parent.field, 'R'));
            self@replab.Rep(parent.group, 'C', parent.dimension, 'isUnitary', parent.isUnitary);
        end

    end

    methods % Implementations

        % Rep

        function s = headerStr(self)
            s = 'Complexification of representation';
        end

        function M = image(self, g, varargin)
            M = self.parent.image(g, varargin{:});
        end

        function M = inverseImage(self, g, varargin)
            M = self.parent.inverseImage(g, varargin{:});
        end

        function b = canComputeType(self, type)
            b = self.parent.canComputeType(type);
        end

    end

end
