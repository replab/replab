classdef YoungDiagram < replab.Str

    properties (SetAccess = protected)
        partition % (`+replab.Partition`): Partition corresponding to the Young diagram
    end

    methods

        function self = YoungDiagram(partition)
            self.partition = partition;
        end

        function [above, left] = aboveLeft(self)
        % Calculates positional information of the Young diagram
        %
        % Returns
        % -------
        % above:
        %   integer(1,\*): Index of the box immediately to the top in the Young diagram, 0 if none
        % left:
        %   integer(1,\*): Index of the box immediately on the left in the Young diagram 0 if none
            part = self.partition;
            n = sum(part);
            m = numel(part);
            cSum = [0 cumsum(part(1:m-1))];
            above = zeros(1, n);
            left = 0:(n-1);
            left(cSum + 1)=0;
            for j = 2:m
                inds = 1:part(j);
                above(cSum(j) + inds) = cSum(j-1) + inds;
            end
        end

    end

end
