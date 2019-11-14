classdef FullySymmetricShape < replab.tensor.Shape
% Describes the shape of a tensor fully symmetric in all its indices
    
    properties (SetAccess = immutable)
        n % integer: tensor rank 
        d % integer: dimension
        dimTable % integer matrix: dimensions of blocks
    end
    
    methods
       
        function self = FullySymmetricShape(n, d)
        % Constructor
        %
        % Args:
        %   n (integer): Tensor rank, number of indices
        %   d (integer): Dimension of each of the indices (which is uniform)
            dimensions = ones(1, n) * d;
            group = replab.Permutations(n);
            self = self@replab.tensor.Shape(dimensions, group, true);
            self.n = n;
            self.d = d;
            dimTable = zeros(d, n);
            for i = 1:d
                for j = 1:n
                    dimTable(i, j) = replab.tensor.FullySymmetricShape.computeDimension(i, j);
                end
            end
            self.dimTable = dimTable;
        end
        
       function p = subOrbit(self, sub)
           [u, ~, J] = unique(sub);
           % Code adapted from https://www.mathworks.com/matlabcentral/newsreader/view_thread/164470
           p = u(up(J, length(sub)));
           function p = up(J, n)
               ktab = histc(J,1:max(J));
               l = n;
               p = zeros(1, n);
               s = 1;
               for i = 1:length(ktab)
                   k = ktab(i);
                   c = nchoosek(1:l, k);
                   m = size(c,1);
                   [t, ~] = find(~p.');
                   t = reshape(t, [], s);
                   c = t(c,:)';
                   s = s*m;
                   r = repmat((1:s)',[1 k]);
                   q = accumarray([r(:) c(:)], i, [s n]);
                   p = repmat(p, [m 1]) + q;
                   l = l - k;
               end
           end
       end

       function sub = indToSub(self, ind)
           ind = ind(:);
           nRows = length(ind);
           switch self.n
             case 0
               sub = zeros(nRows, 0);
             case 1
               sub = ind;
             otherwise
               [sortedInd, index] = sort(ind);
               nRows = length(ind);
               sortedSub = replab.tensor.FullySymmetricShape.indToSubSymHelper(sortedInd, self.d, self.n, self.dimTable);
               sub = zeros(nRows, self.n);
               sub(index, :) = sortedSub;
           end
       end
       
       function ind = subToIndSym(self, sub)
           nRows = size(sub, 1);
           switch self.n
             case 0
               ind = ones(nRows, 1);
             case 1
               ind = sub;
             otherwise
               sub = sort(sub')'; % put each row of subindices in increasing order
                                  % sort treats columns individually, so transpose
               [sortedSub, index] = sortrows(sub); % sorts the subindices for speed
                                                   % we have sortedSub = sub(index, :)
               rStart = 1;
               sortedInd = replab.tensor.FullySymmetricShape.subToIndSymHelper(sortedSub, self.d, self.n, self.dimTable);
               ind = zeros(nRows, 1);
               ind(index) = sortedInd;
           end
       end

   end
   
   methods (Static)
      
       function dim = computeDimension(d, n)
            dim = nchoosek(n + d - 1, n);
        end

        function sub = indToSubSymHelper(ind, d, n, dimTable)
        % helper function for indToSubSym
        % assumes that ind is already sorted
            if n == 1 % handle the trivial case
                sub = ind;
                return
            end
            nRows = length(ind);
            rStart = 1;
            nElementsBefore = 0;
            sub = zeros(nRows, n);
            % for speed, we select the block of indices whose first subindex value
            % is the same; then we call indToSubSymHelper recursively for the remaining
            % columns of the block
            for firstCol = 1:d
                % when the first column index is i, the elements of the remaining
                % columns can be chosen from i to d, thus there are (d - i + 1)
                % choices for these (n - 1) columns
                sizeOfBlock = replab.tensor.FullySymmetricShape.computeDimension(d - firstCol + 1, n - 1);
                startIndex = nElementsBefore + 1;
                endIndex = startIndex + sizeOfBlock - 1;
                if rStart <= nRows && ind(rStart) <= endIndex
                    % we have a block
                    rNextStart = rStart + 1;
                    while rNextStart <= nRows && ind(rNextStart) <= endIndex
                        rNextStart = rNextStart + 1;
                    end
                    rEnd = rNextStart - 1;
                    remainingColsInd = ind(rStart:rEnd) - nElementsBefore;
                    remainingColsSub = replab.tensor.FullySymmetricShape.indToSubSymHelper(remainingColsInd, d - firstCol + 1, n - 1, dimTable);
                    sub(rStart:rEnd, 1) = firstCol;
                    sub(rStart:rEnd, 2:end) = remainingColsSub + firstCol - 1;
                    rStart = rNextStart;
                end
                nElementsBefore = nElementsBefore + sizeOfBlock;
            end            
        end
        
        function ind = subToIndSymHelper(sub, d, n, dimTable)
        % helper function for subToIndSym
        % assumes that sub is already sorted, i.e.
        % each subindex is in the canonical form (increasing)
        % and the rows are sorted lexicographically
            if n == 1
                ind = sub;
                return
            end
            rStart = 1;
            nRows = size(sub, 1);
            ind = zeros(nRows, 1);
            % for speed, we treat all rows with the same "first column value" as a group
            % then, the subspace spanned by the remaining column is also a symmetric subspace
            % whose dimension is reduced by the index of the first column value (canonical indices
            % are increasing), and number of copies is n - 1.
            while rStart <= nRows
                rNextStart = rStart + 1;
                firstCol = sub(rStart, 1);
                while rNextStart <= nRows && sub(rNextStart, 1) == firstCol
                    rNextStart = rNextStart + 1;
                end
                rEnd = rNextStart - 1;
                nElementsBefore = 0;
                for firstColBefore = 1:(firstCol-1)
                    % when the first column index is i, the elements of the remaining
                    % columns can be chosen from i to d, thus there are (d - i + 1)
                    % choices for these (n - 1) columns
                    sizeOfBlock = dimTable(d - firstColBefore + 1, n - 1);
                    nElementsBefore = nElementsBefore + sizeOfBlock;
                end
                remainingCols = sub(rStart:rEnd, 2:end) - firstCol + 1;
                remainingColsInd = replab.tensor.FullySymmetricShape.subToIndSymHelper(remainingCols, d - firstCol + 1, n - 1, dimTable);
                ind(rStart:rEnd) = nElementsBefore + remainingColsInd;
                rStart = rNextStart;
            end
        end

   end
   
end
