classdef Monoid < handle
% Describes a monoid
    
    properties (SetAccess = protected)
        identity % Monoid identity element
    end
    
    methods % Abstract methods
        
        function z = compose(self, x, y)
        % Composes two monoid/group elements
        %
        % Args:
        %   x (element): Left hand side of the binary operation
        %   y (element): Right hand side of the binary operation
        %
        % Returns:
        %   element: Result of the binary operation ``x`` op ``y``
            error('Abstract');
        end

    end
    
end
