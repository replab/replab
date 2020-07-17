classdef PermutationTo1 < replab.bsgs.Backtrack1
% Computes the permutation that sends a vector to another vector
%
% We return the set of ``p`` such that ``t == s(inverse(p))`` or ``s == t(p)``.
%

    properties
        partition % (`+replab.Partition`): Partition to stabilize
        blockIndex % (integer(1,domainSize)): Block index for each point, non-decreasing along the base
    end

    methods

        function self = OrderedPartitionStabilizer1(group, partition, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = group.trivialSubgroup;
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
            base = [blocks{:}];
            self@replab.bsgs.Backtrack1(group.chain, base, knownSubgroup.chain, debug);
            self.partition = partition;
            self.blockIndex = partition.blockIndex;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = self.blockIndex(beta) == self.blockIndex(b);
        end


        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = all(self.blockIndex == self.blockIndex(g));
        end

    end

end
