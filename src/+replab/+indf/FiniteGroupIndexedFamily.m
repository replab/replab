classdef FiniteGroupIndexedFamily < replab.IndexedFamily
% Describes an indexed family of permutations, where the permutations are stored explicitly in a matrix
%
% The permutations are hashed with a simple scheme to faciliate fast retrieval.

    properties (SetAccess = protected)
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        permSet % (`+replab.+perm.Set`): Underlying set of permutations
    end

    methods

        function self = FiniteGroupIndexedFamily(matrix, isomorphism)
        % Constructs an indexed family of finite group elements
        %
        % Args:
        %   matrix (domainSize,\*): Matrix of permutation realizations
        %   isomorphism (`+replab.FiniteIsomorphism`, optional): Isomorphism to a permutation group
            ds = size(matrix, 1);
            if nargin < 2
                isomorphism = [];
            end
            ps = replab.perm.Set(ds);
            ps.insert(matrix);
            self.permSet = ps;
            self.isomorphism = isomorphism;
            self.nElements = vpi(size(matrix, 2));
        end

        function obj = at(self, ind)
            obj = self.permSet.at(double(ind))';
            if ~isempty(self.isomorphism)
                obj = self.isomorphism.preimageElement(obj);
            end
        end

        function ind = find(self, obj)
            if ~isempty(self.isomorphism)
                obj = self.isomorphism.imageElement(obj);
            end
            ind = vpi(self.permSet.find(obj'));
        end

        function C = toCell(self)
            matrix = self.permSet.matrix;
            isomorphism = self.isomorphism;
            if isempty(self.isomorphism)
                C = arrayfun(@(i) matrix(:,i)', 1:size(matrix,2), 'uniform', 0);
            else
                C = arrayfun(@(i) self.isomorphism.preimageElement(matrix(:,i)'), 1:size(matrix,2), 'uniform', 0);
            end
        end

    end

end
