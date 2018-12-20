classdef Partition < replab.Str
   
    properties (SetAccess = protected)
        n;           % domain size
        blockIndex;  % index of the block for each element
        start;       % starting index for each block
        next;        % next index in the same block, or 0 if at the end
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
        
        function s = str(self)
            s = 'Partition ';
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
        
        function B = blocks(self)
            nB = self.nBlocks;
            B = cell(1, nB);
            for i = 1:nB
                B{i} = self.block(i);
            end
        end
        
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
                        c = (blockIndex == 0) & (adjacencyMatrix(t, :) | adjacencyMatrix(:, t)');
                        newAdded = find(c);
                        added = [added newAdded];
                        blockIndex(newAdded) = block;
                        test = [test newAdded];
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
