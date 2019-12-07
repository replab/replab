classdef Partition < replab.Str
% Represents an unordered partition of the set {1..n} into disjoint subsets
   
    properties (SetAccess = protected)
        n % integer: Domain size
        blockIndex % integer row vector: Index of the block for each element
        start % integer row vector: Starting index for each block
        next % integer row vector: Next index in the same block, or 0 if at the end
    end
    
    methods (Access = protected)
        
        function self = Partition(n, blockIndex, start, next)
            self.n = n;
            self.blockIndex = blockIndex;
            self.start = start;
            self.next = next;
        end
        
    end
    
    methods
        
        function s = shortStr(self, maxColumns)
            s = '';
            for i = 1:self.nBlocks
                if i > 1
                    s = [s '|'];
                end
                b = self.block(i);
                for j = 1:length(b)
                    if j > 1 && self.n > 9
                        s = sprintf('%s %d', s, b(j));
                    else
                        s = sprintf('%s%d', s, b(j));
                    end
                end
            end
        end
        
        function lines = longStr(self, maxRows, maxColumns)
            lines = replab.str.longStr(self, maxRows, maxColumns);
            lines{1} = ['Partition ''' self.shortStr(maxColumns) ''''];
        end
        
        function n = nBlocks(self)
            n = length(self.start);
        end
        
        function B = block(self, i)
            el = self.start(i);
            B = [];
            while el > 0
                B = [B el];
                el = self.next(el);
            end
        end
        
        function sz = blockSizes(self)
        % Returns the sizes of all blocks
        %
        % Returns:
        %   (row integer vector): block sizes
            nB = self.nBlocks;
            sz = arrayfun(@(i) self.blockSize(i), 1:nB);
        end
            
        function sz = blockSize(self, i)
        % Returns the size of a partition block
        %
        % Args:
        %   i (integer): Index of the block
        %
        % Returns:
        %   integer: Size of the i-th block in this partition
            el = self.start(i);
            sz = 0;
            while el > 0
                sz = sz + 1;
                el = self.next(el);
            end
        end
        
        function B = blocks(self)
            nB = self.nBlocks;
            B = cell(1, nB);
            for i = 1:nB
                B{i} = self.block(i);
            end
        end
        
        function [P1 pind] = restrictedToBlocks(self, blocks)
        % Returns the partition containing only the given blocks,
        % where the selected blocks are ordered
            pind = [];
            n1 = 0;
            blockIndex1 = [];
            start1 = [];
            next1 = [];
            b1 = 1;
            for b = blocks
                block = self.block(b);
                m = length(block);
                blockIndex1 = [blockIndex1 b1 * ones(1, m)];
                start1 = [start1 (n1 + 1)];
                next1 = [next1 (n1+(2:m)) 0];
                pind = [pind block];
                b1 = b1 + 1;
                n1 = n1 + m;
            end
            rest = setdiff(1:self.n, pind);
            pind = [pind rest];
            P1 = replab.Partition(n1, blockIndex1, start1, next1);
        end

        
% $$$         function [P1 blockIndices p] = subPartitionForBlockMask(self, blockMask)
% $$$             blocks = find(blockMask);
% $$$             blockIndices = find(ismember(self.blockIndex, blocks));
% $$$             rest = setdiff(1:self.n, blockIndices);
% $$$             p = [blockIndices rest]; % original from sub
% $$$             pI(p) = 1:self.n; % sub from original
% $$$             pb = [blocks setdiff(1:self.nBlocks, blocks)]; % original block from sub
% $$$             pbI(pb) = 1:self.nBlocks; % sub block from original
% $$$             n1 = length(blockIndices);
% $$$             blockIndex1 = pbI(self.blockIndex(blockIndices));
% $$$             start1 = pI(self.start(blocks));
% $$$             next1 = zeros(1, n1);
% $$$             mask = self.next(blockIndices) > 0;
% $$$             next1(mask) = pI(self.next(p(mask)));
% $$$             P1 = replab.Partition(n1, blockIndex1, start1, next1);
% $$$         end
% $$$         
% $$$         function P1 = permutationLeftAction(self, g)
% $$$         % Permutes the indices of this permutation
% $$$             gI(g) = 1:self.n;
% $$$             blockIndex1 = self.blockIndex(gI);
% $$$             
% $$$             
% $$$         end
% $$$         
% $$$         function [P1 perm] = subPartitionForIndices(self, indices)
% $$$         % perm: indexIntoSubPartition -> indexIntoOriginalPartition
% $$$             blocks = unique(blockIndex(indices));
% $$$             [P1 blockIndices p] = self.subPartitionForBlocks(blocks);
% $$$             %                               p: ^blockIndices -> ^original
% $$$             [s1, p1] = sort(indices);
% $$$             [s2, p2] = sort(blockIndices);
% $$$             assert(isequal(s1, s2), 'The given indices do not match blocks');
% $$$             % indices(p1) is sorted      - p1: sorted -> ^indices
% $$$             % blockIndices(p2) is sorted - p2: sorted -> ^blockIndices
% $$$             pI1(p1) = 1:length(indices); %     ^indices -> sorted
% $$$             pI2(p2) = 1:length(indices); %     ^blockIndices -> sorted
% $$$             % 
% $$$             
% $$$         end
        
    end

    methods (Static)

        function P = fromBlockIndices(blockIndex)
            n = length(blockIndex);
            nBlocks = max(blockIndex);
            start = zeros(1, nBlocks);
            next = zeros(1, n);
            for i = 1:nBlocks
                blockInd = find(blockIndex == i);
                assert(length(blockInd) > 0, 'Blocks cannot be empty');
                start(i) = blockInd(1);
                for j = 1:length(blockInd)-1
                    next(blockInd(j)) = blockInd(j+1);
                end
            end
            P = replab.Partition(n, blockIndex, start, next);
        end
                            
        function P = connectedComponents(adjacencyMatrix)
        % Given an adjacency matrix adj, returns the sets of
        % vertices corresponding to connected components
        %
        % For adj = [0 0 1; 0 0 0; 1 0 0], it returns the partition {[1 3] [2]}
        
            n = size(adjacencyMatrix, 1);
            assert(size(adjacencyMatrix, 2) == n);
            
            edges = replab.graph.adj2edge(adjacencyMatrix);
            [subsets blockIndex start next] = replab.graph.connectedComponents(edges);
            
            % We don't want sparse objects here
            blockIndex = full(blockIndex);
            next = full(next);
            
            % If some elements are isolated, we add them
            connectedVertices = [subsets{:}];
            isolatedVertices = setdiff(1:n, connectedVertices);
            nbConnectedSets = length(subsets);
            
            if length(isolatedVertices) >= 1
                % allocate memory
                subsets{nbConnectedSets + length(isolatedVertices)} = 0;
                start(nbConnectedSets + length(isolatedVertices)) = start(end);
                next(n) = next(end);
                
                % assign values
                co = nbConnectedSets;
                for i = 1:length(isolatedVertices)
                    co = co + 1;
                    subsets{co} = isolatedVertices(i);
                    blockIndex(isolatedVertices(i)) = co;
                    start(co) = isolatedVertices(i);
                end
            end
            
            % construct the Partition object
            P = replab.Partition(n, blockIndex, start, next);
        end
        
        function P = permutationsOrbits(permutations)
        % Returns the partition of the domain 1...N into orbits
        % where permutations are a nG x domainSize double matrix
            n = size(permutations, 2);
            nG = size(permutations, 1);
            blockIndex = zeros(1, n);
            start = [];
            next = zeros(1, n);
            block = 1;
            for i = 1:n
                if blockIndex(i) == 0
                    % new block discovered, starting with i
                    start = [start i];
                    added = i;
                    blockIndex(i) = block;
                    test = i;
                    while ~isempty(test)
                        t = test(1);
                        test = test(2:end);
                        for j = 1:nG
                            ti = permutations(j, t);
                            if blockIndex(ti) == 0
                                added = [added ti];
                                blockIndex(ti) = block;
                                test = [test ti];
                            end
                        end
                    end
                    added = sort(added);
                    for i = 1:length(added) - 1
                        next(added(i)) = added(i + 1);
                    end
                    block = block + 1;
                end
            end
            P = replab.Partition(n, blockIndex, start, next);
        end
                
    end

end
