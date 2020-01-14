classdef Partition < replab.Str
% Represents an unordered partition of the set {1..n} into disjoint subsets

    properties (SetAccess = protected)
        n % integer: Domain size
        blockIndex % integer row vector: Index of the block for each element
        start % integer row vector: Starting index for each block
        next % integer row vector: Next index in the same block, or 0 if at the end
        blocks % cell array row vector of integer row vector: group elements by partition
    end

    methods (Access = protected)

        function self = Partition(n, blockIndex, start, next, blocks)
            self.n = n;
            self.blockIndex = blockIndex;
            self.start = start;
            self.next = next;
            self.blocks = blocks;
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
            B = self.blocks{i};
        end

        function sz = blockSizes(self)
        % Returns the sizes of all blocks
        %
        % Returns:
        %   (row integer vector): block sizes
            nB = self.nBlocks;
            sz = arrayfun(@(i) length(self.blocks{i}), 1:nB);
        end

        function sz = blockSize(self, i)
        % Returns the size of a partition block
        %
        % Args:
        %   i (integer): Index of the block
        %
        % Returns:
        %   integer: Size of the i-th block in this partition
            sz = length(self.block{i});
        end

        function [P1 pind] = restrictedToBlocks(self, selBlocks)
        % Returns the partition containing only the given blocks
        %
        % The selected blocks are ordered
            pind = [];
            n1 = 0;
            blockIndex1 = [];
            start1 = [];
            next1 = [];
            b1 = 1;
            newBlocks = cell(1,length(selBlocks));
            for b = selBlocks
                block = self.block(b);
                m = length(block);
                blockIndex1 = [blockIndex1 b1 * ones(1, m)];
                start1 = [start1 (n1 + 1)];
                next1 = [next1 (n1+(2:m)) 0];
                pind = [pind block];
                newBlocks{b1} = n1+[1:m];
                b1 = b1 + 1;
                n1 = n1 + m;
            end
            rest = setdiff(1:self.n, pind);
            pind = [pind rest];
            P1 = replab.Partition(n1, blockIndex1, start1, next1, newBlocks);
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

            % Construct the subsets
            blocks = cell(1, nBlocks);
            for i = 1:nBlocks
                blocks{i} = find(blockIndex == i);
                assert(length(blocks{i}) > 0, 'Blocks cannot be empty');
            end

            % Construct the start vector
            start = zeros(1,length(blocks));
            for i = 1:length(blocks)
                start(i) = blocks{i}(1);
            end

            % Construct the next vector
            c = zeros(1,length(a)-length(blocks));
            d = zeros(size(c));
            co = 0;
            for i = 1:length(blocks)
                for j = 1:length(blocks{i})-1
                    co = co + 1;
                    c(co) = blocks{i}(j);
                    d(co) = blocks{i}(j+1);
                end
            end
            next = full(sparse(1,c,d,1,max(max(edges))));

            % Construct the Partition object
            P = replab.Partition(n, blockIndex, start, next, blocks);
        end

        function P = connectedComponentsFromEdges(edges, n)
        % Given list of edges, returns the sets of vertices corresponding to connected components
        %
        % For edges = [1 3] and n = 3, it returns the partition {[1 3] [2]}

            if isempty(edges)
                % Trivial case
                blockIndex = 1:n;
                start = 1:n;
                next = zeros(1,n);
                blocks = num2cell(1:n, 1);
            else
                assert(max(edges(:)) <= n);
                assert(size(edges,2) == 2);

                [blocks blockIndex start next] = replab.graph.connectedComponents(edges);

                % We don't want sparse objects here
                blockIndex = full(blockIndex);
                next = full(next);

                % If some elements are isolated, we add them
                connectedVertices = [blocks{:}];
                isolatedVertices = setdiff(1:n, connectedVertices);
                nbConnectedSets = length(blocks);

                if length(isolatedVertices) >= 1
                    % allocate memory
                    blocks{nbConnectedSets + length(isolatedVertices)} = 0;
                    start(nbConnectedSets + length(isolatedVertices)) = start(end);
                    next(n) = next(end);

                    % assign values
                    co = nbConnectedSets;
                    for i = 1:length(isolatedVertices)
                        co = co + 1;
                        blocks{co} = isolatedVertices(i);
                        blockIndex(isolatedVertices(i)) = co;
                        start(co) = isolatedVertices(i);
                    end
                end
            end

            % Construct the Partition object
            P = replab.Partition(n, blockIndex, start, next, blocks);
        end

        function P = connectedComponents(adjacencyMatrix)
        % Given an adjacency matrix adj, returns the sets of vertices corresponding to connected components
        %
        % For adj = [0 0 1; 0 0 0; 1 0 0], it returns the partition {[1 3] [2]}

            n = size(adjacencyMatrix, 1);
            assert(size(adjacencyMatrix, 2) == n);

            edges = replab.graph.adj2edge(adjacencyMatrix);
            P = replab.Partition.connectedComponentsFromEdges(edges, n);
        end

        function P = permutationsOrbits(permutations)
        % Returns the partition of the domain 1...N into orbits
        %
        % The permutations are a nG x domainSize double matrix
            n = size(permutations, 2);
            nG = size(permutations, 1);

            % We list the edges of the graph
            edges = cell(nG,1);
            for i = 1:nG
                edges{i} = [(1:n); permutations(i,:)].';
            end
            edges = unique(cat(1, edges{:}), 'rows');

            % Call connected component method to construct the Partition
            % object
            P = replab.Partition.connectedComponentsFromEdges(edges, n);
        end

    end

end
