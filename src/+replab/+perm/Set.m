classdef Set < replab.Str
% Stores a set of permutations
%
% This is not a multiset; permutations are stored only once, and the implementation is mutable.
%
% Those permutations are stored, not necessarily ordered, as columns in a matrix.
% In parallel, we compute a hash of each of those permutations, and we stored
% a sorted list of hash values, and the relation between the ordering of permutations
% in this sorted and the ordering in the matrix of permutations.
%
% This enables the use of fast vectorized operations in Matlab while retrieving elements.

    properties
        domainSize % (integer): Size of the stored permutations
        matrix % (integer(domainSize,size)): Elements stored as column vectors
        seed % integer(1,domainSize): Vector used to hash a permutation
        perm % integer(1,size): Permutation of the columns of `.matrix`
        hash % integer(1,size): Sorted hash vector ``= seed * matrix(:, perm)``
    end

    methods

        function self = Set(domainSize)
            self.domainSize = domainSize;
            self.matrix = zeros(domainSize, 0);
            % we use a simple hash function that maps permutations to doubles
            % we want that ``seed * perm``, is an integer that fits in a double
            % if the seed coefficients have range in ``[-r, r]``, then the resulting hash
            % is in ``[-r*domainSize^2, r*domainSize^2]``
            r = floor(2^52/domainSize/domainSize) - 1;
            self.seed = randi([-r r], 1, domainSize);
            self.perm = zeros(1, 0);
            self.hash = zeros(1, 0);
        end

        function check(self)
        % Checks the consistency of internal data structures
            assert(isequal(self.hash, self.seed * self.matrix(:, self.perm)));
        end

        function s = nElements(self)
        % Returns the number of elements in this set
            s = size(self.matrix, 2);
        end

        function permutations = at(self, inds)
        % Returns the permutations at the given indices
        %
        % Args:
        %   inds (integer(1,\*)): Indices to retrieve
        %
        % Returns:
        %   integer(domainSize,\*): Permutations as column vectors
            permutations = self.matrix(:, inds);
        end

        function sort(self)
        % Sorts the set
            [matrix1, ind] = sortrows(self.matrix');
            matrix1 = matrix1';
            % matrix1 = matrix(:,ind)
            % hash = seed * matrix(:, perm)
            % hash = seed * matrix1(:, perm1)
            % hash = seed * matrix(:, ind(perm1))
            % thus perm = ind(perm1) which is perm = ind * perm1, perm1 = indInv * perm
            % perm1 = indInv(perm)
            n = length(ind);
            indInv = zeros(1, n);
            indInv(ind) = 1:n;
            perm1 = indInv(self.perm);
            self.matrix = matrix1;
            self.perm = perm1;
            % hash doesn't change
        end

        function inds = insert(self, permutations)
        % Inserts the given permutations in the set
        %
        % Modifies the set in place.
        %
        % Args:
        %   permutations (integer(domainSize,\*)): New permutations to insert, must not be already present
        %
        % Returns:
        %   integer(1,\*): Indices of the inserted permutations
            unhash = zeros(1, self.nElements); % unsorted hash
            unhash(self.perm) = self.hash;
            addHash = self.seed*permutations;
            unhash1 = [unhash addHash]; % new unsorted hash
            matrix1 = [self.matrix permutations];
            assert(isequal(unhash1, self.seed * matrix1));
            [hash1, perm1] = sort(unhash1);
            oldSize = self.nElements;
            self.matrix = matrix1;
            self.perm = perm1;
            self.hash = hash1;
            inds = oldSize+1:self.nElements;
        end

        function ind = findElement(self, col, range)
        % Returns the index of the given column, or 0 if not found
        %
        % Args:
        %   col (integer(\*,1)): Integer vector
        %   range (integer(1,\*), optional): Range in which to search for the element
        %
        % Returns:
        %   integer: Column index, or ``0`` if not found
            if isempty(self.matrix)
                ind = 0;
                return
            end
            if nargin < 3 || isempty(range)
                ind = find(self.matrix(1,:) == col(1));
            else
                ind = range(self.matrix(1, range) == col(1));
            end
            nr = size(self.matrix, 1);
            for r = 2:nr
                switch length(ind)
                  case 0
                    return
                  case 1
                    if any(self.matrix(r:end, ind) ~= col(r:end))
                        ind = 0;
                    end
                    return
                  case 2
                    if all(self.matrix(r:end, ind(1)) == col(r:end))
                        ind = ind(1);
                    elseif all(self.matrix(r:end, ind(2)) == col(r:end))
                        ind = ind(2);
                    else
                        ind = 0;
                    end
                    return
                  otherwise
                    ind = range(self.matrix(r, ind) == col(r));
                end
            end
        end


        function inds = find(self, permutations)
        % Returns the indices of the given permutations, or 0 if not found
        %
        % Args:
        %   permutations(domainSize,\*): Permutations to look for
        %
        % Returns:
        %   integer(1,\*): Indices of the given permutations, or ``0`` where not found
            n = size(permutations, 2);
            inds = zeros(1, n);
            hash = self.hash;
            hs = self.seed * permutations;
            fast = exist('ismembc2') > 0;
            if fast
                % perform binary search; this function returns the last element which matches
                lastInds = ismembc2(hs, hash);
            end
            for i = 1:n
                h = hs(i);
                if fast
                    last = lastInds(i);
                    if last == 0
                        inds(i) = 0;
                        continue
                    end
                    % then we need to find if other elements before match as well
                    before = last - 1;
                    while before > 0 && hash(before) == h
                        before = before - 1;
                    end
                    f = before+1:last;
                else
                    f = find(hash == h);
                    if isempty(f)
                        inds(i) = 0;
                        continue
                    end
                end
                ind = self.perm(f);
                if isempty(ind)
                    inds(i) = 0;
                else
                    inds(i) = self.findElement(permutations(:,i), ind);
                end
            end
        end

        function inds = update(self, permutations)
        % Returns the indices of the given permutations, adding them to the set if necessary
        %
        % Modifies the set in place.
        %
        % Args:
        %   permutations(domainSize,\*): Permutations to look for
        %
        % Returns:
        %   integer(1,\*): Indices of the given permutations
            inds = self.find(permutations);
            mask = (inds == 0);
            inds1 = self.insert(permutations(:, find(mask)));
            inds(mask) = inds1;
        end

    end

end
