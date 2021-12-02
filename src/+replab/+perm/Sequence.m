classdef Sequence < replab.Sequence
% Describes a sequence of permutations, where the permutations are stored explicitly in a matrix

    properties (SetAccess = protected)
        permSet % (`+replab.+perm.Set`): Underlying set of permutations
    end

    methods

        function self = Sequence(matrix)
        % Constructs a sequence of permutations
        %
        % Args:
        %   matrix (domainSize,\*): Matrix of permutations
            ds = size(matrix, 1);
            ps = replab.perm.Set(ds);
            ps.insert(matrix);
            self@replab.Sequence(size(matrix, 2));
            self.permSet = ps;
        end

        function obj = at(self, ind)
            if ischar(ind)
                ind = str2num(ind);
            end
            obj = self.permSet.at(double(ind))';
        end

        function ind = find(self, obj)
            ind = vpi(self.permSet.find(obj'));
        end

        function C = toCell(self)
            matrix = self.permSet.matrix;
            C = arrayfun(@(i) matrix(:,i)', 1:size(matrix,2), 'uniform', 0);
        end

    end

end
