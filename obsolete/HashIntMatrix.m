classdef HashIntMatrix
% An int32 matrix who supports fast lookup of columns
    
    properties
        M;
        mask;
        sortedHashs;
        hashIndices;
    end
    
    methods (Static)
        
        function left = leftBinarySearch(vec, target, right)
        % Finds the leftmost target index
            left = 1;
            if nargin < 3
                right = length(vec);
            end
            % carry right to the greatest number which is less than target
            while left <= right
                mid = floor((left + right)/2);
                if vec(mid) < target
                    left = mid + 1;
                else
                    right = mid - 1;
                end
            end
            % when we are here, right is a the index of the greatest number
            % which is less than target and since left is at the next, it is
            % at the first target's index
        end
        
        function right = rightBinarySearch(vec, target, left)
        % Finds the rightmost target index
            if nargin < 3
                left = 1;
            end
            right = length(vec);
            % carry left to the smallest number which is greater than target
            while left <= right
                mid = floor((left + right)/2);
                if vec(mid) <= target
                    left = mid + 1;
                else
                    right = mid - 1;
                end
            end
            % when we are here, left is at the index of smallest number
            % which is greater than target and since right is at the next, 
            % it is at the first target's index
        end
        
    end
    
    methods
        
        function self = HashIntMatrix(M, mask, sortedHashs, hashIndices)
            M = int32(M);
            if nargin < 4
                d = size(M, 1);
                n = size(M, 2);
                mask = randi(20, 1, d, 'int32')*2 + 1;
                hashs = zeros(1, n, 'int32');
                for i = 1:d
                    hashs = hashs + mask(i) * M(i,:);
                end
                [sortedHashs hashIndices] = sort(hashs);
            end
            self.M = M;
            self.mask = mask;
            self.sortedHashs = sortedHashs;
            self.hashIndices = hashIndices;
        end
        
        function n = nElements(self)
            n = size(self.M, 2);
        end
        
        function newHIM = append(self, newCols)
            newCols = int32(newCols);
            nNew = size(newCols, 2);
            d = size(self.M, 1);
            oldPerm = [self.hashIndices length(self.hashIndices)+(1:nNew)];
            appendHashs = zeros(1, nNew, 'int32');
            for i = 1:d
                appendHashs = appendHashs + self.mask(i) * newCols(i, :);
            end
            newHashs = [self.sortedHashs appendHashs];
            [newSortedHashs p] = sort(newHashs);
            newHashIndices = oldPerm(p);
            newM = [self.M newCols];
            newHIM = replab.prv.HashIntMatrix(newM, self.mask, newSortedHashs, newHashIndices);
        end
        
        function [newHIM, indices] = lexColSorted(self)
            [newM, indices] = sortrows(self.M');
            invIndices = replab.Perm.inverse(indices);
            newHIM = replab.prv.HashIntMatrix(newM');
            newM = self.M(:, indices);
            newHashIndices = invIndices(self.hashIndices);
            newHIM = replab.prv.HashIntMatrix(newM, self.mask, self.sortedHashs, newHashIndices);
        end
        
        function ind = find(self, targetCol)
        % Finds the matrix col equal to targetCol
            ind = 0;
            targetCol = int32(targetCol(:));
            d = length(self.mask);
            hash = int32(0);
            for i = 1:d
                hash = hash + self.mask(i) * targetCol(i);
            end
            left = replab.prv.HashIntMatrix.leftBinarySearch(self.sortedHashs, hash);
            if left < 1 || left > self.nElements || self.sortedHashs(left) ~= hash
                return
            end
            right = replab.prv.HashIntMatrix.rightBinarySearch(self.sortedHashs, hash, left);
            for c = self.hashIndices(left:right)
                if isequal(self.M(:,c), targetCol)
                    ind = c;
                    return
                end
            end
        end
        
    end
    
end
