classdef TotalOrder < replab.Domain
% Describes a domain equipped with a total order

    methods % Comparison

        function c = compare(self, x, y)
        % Compares two elements of a domain
        %
        % The method returns:
        %
        % - ``-1`` if ``x < y``
        % - ``0`` if ``x == y``
        % - ``1`` if ``x > y``
        %
        % Args:
        %   x (element): First element
        %   y (element): Second element
        %
        % Returns:
        %   -1,0,1: Result of comparison
            error('Abstract');
        end

    end

end
