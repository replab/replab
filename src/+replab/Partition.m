classdef Partition < replab.Str
% Represents an unordered partition of the set ``{1..n}`` into disjoint subsets
%
% The subsets, or blocks, are represented by sorted integer row vectors. The subsets
% themselves are sorted by their minimal element.

    properties (SetAccess = protected)
        n % (integer): Domain size
        blockIndex % (integer(1, n)): Index of the block for each element
        blocks % (cell(1,\*) of integer(1,\*)): List of blocks
    end

    methods (Access = protected)

        function self = Partition(blockIndex, blocks)
            self.n = length(blockIndex);
            self.blockIndex = blockIndex;
            self.blocks = blocks;
        end

    end

    methods

        function check(self)
        % Verifies the sanity of this partition
            m = cellfun(@min, self.blocks); % blocks are ordered
            assert(all(m(2:end) - m(1:end-1)) > 0);
            % all elements accounted for
            assert(sum(cellfun(@length, self.blocks)) == self.n);
            for i = 1:self.nBlocks
                block = self.block(i);
                assert(isequal(block, sort(block))); % each block is a sorted vector
                assert(self.blockIndex(block) == i); % blockIndex and blocks are consistent
            end
        end

        function l = ne(self, rhs)
        % Checks if this partition differs to another partition
        %
        % Args:
        %   rhs (`.Partition`): Another partition
        %
        % Returns:
        %   logical: True is both partitions differ
            l = ~(self == rhs);
        end

        function l = eq(self, rhs)
        % Checks if this partition is equal to another partition
        %
        % Args:
        %   rhs (`.Partition`): Another partition
        %
        % Returns:
        %   logical: True is both partitions are equal
            l = isequal(self.blockIndex, rhs.blockIndex);
        end

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
        % Returns the number of subsets/blocks in this partition
        %
        % Returns:
        %   integer: Number of blocks
            n = length(self.blocks);
        end

        function B = block(self, i)
        % Returns the ``i``-th subset in this partition
        %
        % Args:
        %   i (integer): Block index
        %
        % Returns:
        %   integer(1,\*): Subset
            B = self.blocks{i};
        end

        function sz = blockSizes(self)
        % Returns the sizes of blocks
        %
        % Returns:
        %   integer(1,\*): block sizes
            nB = self.nBlocks;
            sz = arrayfun(@(i) length(self.blocks{i}), 1:nB);
        end

        function sz = blockSize(self, i)
        % Returns the size of the ``i``-th block
        %
        % Args:
        %   i (integer): Block index
        %
        % Returns:
        %   integer: Size of the ``i``-th block in this partition
            sz = length(self.block{i});
        end

        function [P1 pind] = restrictedToBlocks(self, selBlocks)
        % Returns the partition containing only the given blocks
        %
        % The selected blocks are ordered
            pind = [];
            n1 = 0;
            blockIndex1 = [];
            b1 = 1;
            blocks1 = cell(1,length(selBlocks));
            for b = selBlocks
                block = self.block(b);
                m = length(block);
                blockIndex1 = [blockIndex1 b1 * ones(1, m)];
                pind = [pind block];
                blocks1{b1} = n1+[1:m];
                b1 = b1 + 1;
                n1 = n1 + m;
            end
            rest = setdiff(1:self.n, pind);
            pind = [pind rest];
            P1 = replab.Partition(blockIndex1, blocks1);
        end

    end

    methods (Static)

        function P = fromBlocks(blocks)
        % Constructs a partition from disjoint blocks
        %
        % Example:
        %   >>> replab.Partition.fromBlocks({[1 2 5] [3 4]})
        %     Partition '125|34'
        %     blockIndex: [1, 1, 2, 2, 1]
        %         blocks: {[1, 2, 5], [3, 4]}
        %              n: 5
        %
        % Args:
        %   blocks (cell(1,\*) of integer(1,\*)): Disjoint blocks
        %
        % Returns:
        %   `+replab.Partition`: Constructed partition
            blocks = cellfun(@(b) sort(b), blocks, 'uniform', 0);
            numEl = cellfun(@(b) length(b), blocks);
            minEl = cellfun(@(b) min(b), blocks);
            maxEl = cellfun(@(b) max(b), blocks);
            [~, I] = sort(minEl);
            blocks = blocks(I);
            n = sum(numEl);
            assert(n == max(maxEl));
            assert(1 == min(minEl));
            blockIndex = zeros(1, n);
            for i = 1:length(blocks)
                blockIndex(blocks{i}) = i;
            end
            P = replab.Partition(blockIndex, blocks);
        end

        function P = fromVector(vec)
        % Returns the partition that groups equal coefficients of a vector
        %
        % Example:
        %   >>> replab.Partition.fromVector([0 0 1 1 0])
        %     Partition '125|34'
        %     blockIndex: [1, 1, 2, 2, 1]
        %         blocks: {[1, 2, 5], [3, 4]}
        %              n: 5
        %
        % Args:
        %   vec (double(1,\*)): Vector to group the coefficients of
        %
        % Returns:
        %   `.replab.Partition`: Partition of blocks with equal coefficients
            assert(isvector(vec));
            v = unique(vec);
            blocks = {};
            for i = 1:length(v)
                blocks{1, end+1} = find(vec == v(i));
            end
            P = replab.Partition.fromBlocks(blocks);
        end

        function P = connectedComponentsFromEdges(edges, n)
        % Given list of edges, returns the sets of vertices corresponding to connected components
        %
        % For edges = [1 3] and n = 3, it returns the partition {[1 3] [2]}

            if isempty(edges)
                % Trivial case
                blockIndex = 1:n;
                blocks = num2cell(1:n, 1);
            else
                assert(max(edges(:)) <= n);
                assert(size(edges,2) == 2);

                [blocks blockIndex] = replab.graph.connectedComponents(edges);

                % If some elements are isolated, we add them
                connectedVertices = [blocks{:}];
                isolatedVertices = setdiff(1:n, connectedVertices);
                nbConnectedSets = length(blocks);

                if length(isolatedVertices) >= 1
                    % allocate memory
                    blocks{nbConnectedSets + length(isolatedVertices)} = [];
                    blockIndex(n) = 0;

                    % assign values
                    co = nbConnectedSets;
                    for i = 1:length(isolatedVertices)
                        co = co + 1;
                        blocks{co} = isolatedVertices(i);
                        blockIndex(isolatedVertices(i)) = co;
                    end
                end
            end

            % Construct the Partition object
            P = replab.Partition(blockIndex, blocks);
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
        % Returns the partition of the domain ``1...N`` into orbits
        %
        % Args:
        %   permutations (integer(nG, d)): Permutations given as rows in a matrix
            d = size(permutations, 2);
            nG = size(permutations, 1);
            blockIndex = zeros(1, d);
            nBlocks = 1;
            i = 1;
            while ~isempty(i)
                toCheck = i;
                blockIndex(toCheck) = nBlocks;
                while ~isempty(toCheck)
                    images = permutations(:, toCheck);
                    images = images(:);
                    toCheck = images(blockIndex(images) == 0);
                    toCheck = toCheck(:).';
                    blockIndex(toCheck) = nBlocks;
                end
                nBlocks = nBlocks + 1;
                i = find(blockIndex(i+1:end) == 0, 1) + i;
            end
            P = replab.Partition.fromVector(blockIndex);
        end

    end

end
