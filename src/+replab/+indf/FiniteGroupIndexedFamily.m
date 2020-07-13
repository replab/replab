classdef FiniteGroupIndexedFamily < replab.IndexedFamily
% Describes an indexed family of permutations, where the permutations are stored explicitly in a matrix
%
% The permutations are hashed with a simple scheme to faciliate fast retrieval.

    properties (SetAccess = protected)
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        matrix % (integer(\*,\*)): Elements stored as column vectors
        seed % integer(1, \*): Vector used to hash a permutation
        perm % integer(1, \*): Permutation of the columns of `.matrix`
        hash % integer(1, \*): Sorted hash vector ``= seed * matrix(:, perm)``
    end

    methods

        function self = FiniteGroupIndexedFamily(matrix, isomorphism)
            ds = size(matrix, 1);

            if nargin < 2
                isomorphism = [];
            end
            self.size = vpi(size(matrix, 2));
            self.matrix = matrix;
            if ~isempty(matrix)
                ds = size(matrix, 1);
                n = size(matrix, 2);
                r = floor(2^52/n/n) - 1;
                % we use a simple hash function that maps permutations to doubles
                % with a domain size < 2^16, that means that if h has values between -2^20+1 and 2^20-1,
                % the maximal value of the hash is 2^16*(2^16*2^20) = 2^52 which fits in a double
                seed = randi([-r r], 1, ds)-1;
                unsortedHash = seed * matrix;
                [hash, perm] = sort(unsortedHash);
                invPerm = zeros(1, n);
                invPerm(perm) = 1:n;
                self.seed = seed;
                self.perm = perm;
                self.hash = hash;
            end
        end

        function obj = at(self, ind)
            obj = self.matrix(:, ind)';
            if ~isempty(self.isomorphism)
                obj = self.isomorphism.preimageElement(obj);
            end
        end

        function ind = find(self, obj)
            if isempty(self.matrix)
                ind = vpi(0);
                return
            end
            if ~isempty(self.isomorphism)
                obj = self.isomorphism.imageElement(obj);
            end
            hash = self.hash;
            h = self.seed * obj(:);
            if exist('ismembc2') > 0;
                % perform binary search; this function returns the last element which matches
                last = ismembc2(h, hash);
                % then we need to find if other elements before match as well
                before = last - 1;
                while before > 0 && hash(before) == h
                    before = before - 1;
                end
                f = before+1:last;
            else
                f = find(hash == h);
            end
            ind = self.perm(f);
            if length(ind) > 1
                % several rows have the same hash, so we look for an exact match
                loc = replab.util.findRowInMatrix(obj, self.matrix(:,ind)');
                ind = ind(loc);
            end
            if length(ind) == 1
                if ~isequal(obj, self.matrix(:,ind)')
                    ind = vpi(0);
                end
            end
        end

        function C = toCell(self)
            matrix = self.matrix;
            isomorphism = self.isomorphism;
            if isempty(self.isomorphism)
                C = arrayfun(@(i) self.matrix(:,i)', 1:size(matrix,2), 'uniform', 0);
            else
                C = arrayfun(@(i) self.isomorphism.preimageElement(self.matrix(:,i)'), 1:size(matrix,2), 'uniform', 0);
            end
        end

    end

end
