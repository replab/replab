classdef UnorderedPartitionStabilizer1 < replab.bsgs.Backtrack1
% Computes the unordered partition stabilizer of a group

    properties
        partition % (`+replab.Partition`): Partition to stabilize
        sourceBlockIndex % (integer(1,domainSize)): Block index for each point, non-decreasing along the base
        targetBlockIndex % (integer(1,nBlocks)): Target block index for each source block; is initialized when entering the source block
        blocks % (cell(1,nBlocks) of integer(1,\*)): Blocks of the partition, non-decreasing along the base
        blockSize % (integer(1,domainSize)): Block size for each point in the domain
    end

    methods

        function self = UnorderedPartitionStabilizer1(group, partition, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            n = partition.n;
            assert(n == group.domainSize);

            % sort the blocks from the smallest to the biggest
            blocks = partition.blocks;
            lengths = cellfun(@length, blocks);
            mask = lengths > 1; % filter singleton blocks
            blocks = blocks(mask);
            lengths = lengths(mask);
            [~, I] = sort(cellfun(@length, blocks));
            blocks = blocks(I);
            sourceBlockIndex = zeros(1, n);
            base = zeros(1, 0);
            nB = length(blocks);
            blockSize = zeros(1, n);
            for i = 1:nB
                block = blocks{i};
                blockSize(block) = length(block);
                sourceBlockIndex(block) = i;
                base = [base block];
            end

            self@replab.bsgs.Backtrack1(group, base, knownSubgroup, debug);

            self.blocks = blocks;
            self.partition = partition;
            self.sourceBlockIndex = sourceBlockIndex;
            self.targetBlockIndex = zeros(1, nB);
            self.blockSize = blockSize;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            sbi = self.sourceBlockIndex(beta);
            tbi = self.sourceBlockIndex(b); % source and target block index for the current base point
            ok = self.blockSize(beta) == self.blockSize(b); % sizes must match
            if ~ok
                return
            end
            if l == 1 % we are starting a new block
                prevBeta = 0;
                prevSbi = 0;
            else
                prevBeta = self.base(l-1); % we are maybe starting a new block
                prevSbi = self.sourceBlockIndex(prevBeta);
            end
            if prevSbi ~= sbi % new block
                self.targetBlockIndex(sbi) = tbi; % memorize the block index of the image
            else
                ok = self.targetBlockIndex(prevSbi) == tbi; % verify the block index of the image
            end
        end


        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            for i = 1:length(self.blocks)
                block = self.blocks{i};
                tbi = self.sourceBlockIndex(g(block));
                if any(tbi ~= tbi(1))
                    ok = false;
                    return
                end
            end
            ok = true;
        end

    end

end
