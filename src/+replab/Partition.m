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

    methods

        function self = Partition(blockIndex, blocks)
        % Constructs a Partition from blockIndex and list of blocks
        %
        % Do not use this function direclty, rather use another
        % constructor such as ``.fromBlocks``.
        %
        % Args:
        %   blockIndex
        %   blocks (cell(1,\*) of integer(1,\*)): Disjoint blocks
        %
        % Returns:
        %   `.Partition`: The partition
        %
        % See also:
        %   `.fromBlocks`
        %   `.check`

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

    end

    methods % Properties

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

    end

    methods % Binary operations

    end

    methods % Implementations

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
            if ~isa(rhs, 'replab.Partition')
                l = false;
                return
            end
            l = isequal(self.blockIndex, rhs.blockIndex);
        end

        % Str

        function s = shortStr(self, maxColumns)
            s = '';
            for i = 1:min(self.nBlocks, maxColumns)
                if i > 1
                    s = [s '|'];
                end
                b = self.block(i);
                for j = 1:min(length(b), maxColumns)
                    if j > 1 && self.n > 9
                        s = sprintf('%s %d', s, b(j));
                    else
                        s = sprintf('%s%d', s, b(j));
                    end
                end
                if length(b) > maxColumns
                    s = sprintf('%s...', s);
                end
            end
            if self.nBlocks > maxColumns
                s = sprintf('%s...', s);
            end
        end

        function lines = longStr(self, maxRows, maxColumns)
            lines = replab.str.longStr(self, maxRows, maxColumns);
            lines{1} = ['Partition ''' self.shortStr(maxColumns) ''''];
        end

    end

    methods

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

        function S = singletons(self)
        % Returns a set of all points that are singletons of this partition
        %
        % The singletons are the blocks of size 1.
        %
        % Returns:
        %   integer(1,\*): Set of points
            blocks = self.blocks;
            lengths = cellfun(@length, blocks);
            blocks = blocks(lengths == 1);
            S = [blocks{:}];
        end

    end

    methods (Static)

        function P = trivial(n)
        % Constructs the trivial partition with a single block of given cardinality
        %
        % Args:
        %   n (integer): Size of the partition and single block
        %
        % Returns:
        %   `.Partition`: Trivial partition
            P = replab.Partition.fromBlocks({1:n});
        end

        function P = finest(n)
        % Constructs the partition of ``n`` singleton blocks
        %
        % Args:
        %   n (integer): Size of the partition and number of blocks
        %
        % Returns:
        %   `.Partition`: Partition
            P = replab.Partition.fromBlocks(num2cell(1:n));
        end

        function P = fromBlocks(blocks)
        % Constructs a partition from disjoint blocks
        %
        % The blocks will be sorted and reorganized.
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

        function P = fromApproximateVector(vec, tol)
        % Returns the partition that groups approximately equal coefficients of a vector
        %
        % Example:
        %   >>> replab.Partition.fromApproximateVector([0.1 0.2 0.9 1 0.3], 0.11)
        %     Partition '125|34'
        %     blockIndex: [1, 1, 2, 2, 1]
        %         blocks: {[1, 2, 5], [3, 4]}
        %              n: 5
        %
        % Args:
        %   vec (double(1,\*), real): Vector to group the coefficients of
        %   tol (double): Tolerance
        %
        % Returns:
        %   `.Partition`: Partition of blocks with approximately equal coefficients
            assert(isreal(vec));
            [v1, I] = sort(vec);
            close = [(v1(2:end) - v1(1:end-1) <= tol) false];
            next = 1;
            blocks = cell(1, 0);
            while next <= length(vec)
                current = next;
                next = find(~close(current:end), 1) + current; % + 1 - 1
                blocks{1,end+1} = I(current:next-1);
            end
            P = replab.Partition.fromBlocks(blocks);
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
        %   `.Partition`: Partition of blocks with equal coefficients
            assert(isvector(vec));
            v = unique(vec);
            blocks = {};
            for i = 1:length(v)
                blocks{1, end+1} = find(vec == v(i));
            end
            P = replab.Partition.fromBlocks(blocks);
        end

        function P = permutationsOrbits(permutations)
        % Returns the partition of the domain ``1...N`` into orbits
        %
        % Args:
        %   permutations (integer(nG, d)): Permutations given as rows in a matrix
        %
        % Returns:
        %   `.Partition`: Partition containing the orbits
            d = size(permutations, 2);
            nG = size(permutations, 1);

            % We list the edges of the graph
            edges = permutations.';
            edges = [kron(ones(nG,1), (1:d)'), edges(:)];

            % Compute the connected components
            P = replab.UndirectedGraph.fromEdges(edges, d).connectedComponents;
        end

    end

end
