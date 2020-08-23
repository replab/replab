classdef IntegerPartition < replab.Str

    properties (SetAccess = protected)
        n % integer: Number of boxes in the Young diagram
        partition % (integer(1,\*)): Partition corresponding to the Young diagram
    end

    methods (Static, Access = protected)

        function I = enumerate(n, m)
        % Computes all integer partitions of the given number
        %
        % Args:
        %   n (integer): Integer to partition
        %   m (integer): Upper bound on the parts
        %
        % Returns:
        %   cell(1,\*) of integer(1,\*): Partitions
            if nargin < 2 || isempty(m)
                m = n;
            end
            if n == 0
                I = {zeros(1, 0)};
                return
            end
            if m == 1
                I = {ones(1, n)};
                return
            end
            I = cell(1, 0);
            for i = min(m, n):-1:1
                I = horzcat(I, cellfun(@(part) [i part], replab.sym.IntegerPartition.enumerate(n - i, i), 'uniform', 0));
            end
        end

        function [gens orderFactors] = centralizerBlocks(n, start, blockSize, nBlocks)
        % Helper function used to construct the centralizer of a canonical conjugacy class representative
        %
        % This assumes that the representative has ``nBlocks`` of size ``blockSize`` starting at position ``start``.
        %
        % For the partition ``[3 2 2]``, the representative is ``[2 1 4 3 6 7 5]`` corresponding to blocks of size ``[2 2 3]``
        % which is ``[3 2 2]`` flipped. To construct the centralizer of that element, we will call ``centralizerBlocks`` twice,
        % first with ``(start, blockSize, nBlocks)`` equal to ``(1, 2, 2)`` for the ``[2 2]`` part, then equal to ``(5, 3, 1)`` for the
        % last part ``[3]``.
        %
        % Args:
        %   n (integer): Integer being partitioned
        %   start (integer): Start index, first element of the first block
        %   blockSize (integer): Number of elements in a block, corresponds to the part of the integer partition
        %   nBlocks (integer): How many blocks, corresponds to the number of parts in the integer partition that are equal to ``blockSize``
        %
        % Returns
        % -------
        % gens:
        %   cell(1,\*) of permutation: Generators of the subgroup of the centralizer corresponding to the given blocks
        % order:
        %   integer(1,\*): Factors of the order of that subgroup of the centralizer
            orderFactors = [ones(1,nBlocks)*blockSize 2:nBlocks];
            orderFactors = orderFactors(orderFactors ~= 1);
            gens = cell(1, 0);
            blocks = arrayfun(@(i) (start+(i-1)*blockSize):(start+i*blockSize-1), 1:nBlocks, 'uniform', 0);
            if blockSize > 1
                block1 = blocks{1};
                inCycle = 1:n;
                inCycle(block1) = block1([2:blockSize 1]);
                gens{1,end+1} = inCycle;
            end
            if nBlocks > 1
                swap = 1:n;
                swap([blocks{1} blocks{2}]) = [blocks{2} blocks{1}];
                gens{1,end+1} = swap;
            end
            if nBlocks > 2
                cycle = 1:n;
                cycle([blocks{1:end}]) = [blocks{2:end} blocks{1}];
                gens{1,end+1} = cycle;
            end
        end

    end

    methods (Static)

        function I = all(n)
        % Returns all integer partitions for an integer n
        %
        % Args:
        %   n (integer): Integer to partition
        %
        % Returns:
        %   cell(1,\*) of `.IntegerPartition`: Integer partitions
            I = cellfun(@(p) replab.sym.IntegerPartition(p), replab.sym.IntegerPartition.enumerate(n), 'uniform', 0);
        end

    end

    methods

        function self = IntegerPartition(partition)
            self.partition = partition;
            self.n = sum(partition);
        end

        function C = conjugacyClass(self)
        % Returns the conjugacy class of the symmetric group with cycle structure corresponding to the integer partition
            rep = self.minLexPermutation;
            repCent = self.centralizerMinLexPermutation;
            C = replab.ConjugacyClass(replab.S(self.n), rep, repCent);
        end

        function p = minLexPermutation(self)
        % Returns the minimal permutation under lexicographic ordering that has its cycle structure given by this integer partition
        %
        % Returns:
        %   permutation: The corresponding permutation
            part = fliplr(self.partition);
            p = [];
            ind = 1;
            for j = 1:length(part)
                pj = part(j);
                p = [p ind+(1:(pj-1)) ind];
                ind = ind + pj;
            end
        end

        function G = centralizerMinLexPermutation(self)
        % Returns the subgroup of the symmetric group that centralizes `.minLexPermutation`
        %
        % Returns:
        %   `+replab.PermutationGroup`: The corresponding centralizer
            part = fliplr(self.partition);
            j = 1;
            ind = 1;
            gens = cell(1, 0);
            orderFactors = [];
            while j <= length(part) % for each part size
                j1 = find([part(j+1:end) 0] ~= part(j), 1) + j; % find when the next part size starts
                blockSize = part(j);
                nBlocks = j1 - j;
                % computes the wreath product of C(blockSize) by S(nBlocks)
                [gensJ ofJ] = replab.sym.IntegerPartition.centralizerBlocks(self.n, ind, blockSize, nBlocks);
                % append generators and multiply order
                gens = horzcat(gens, gensJ);
                orderFactors = [orderFactors ofJ];
                ind = ind + nBlocks * blockSize;
                j = j1;
            end
            order = replab.util.multiplyIntegers(orderFactors);
            % construct the group
            G = replab.PermutationGroup(self.n, gens, order, replab.SymmetricGroup.make(self.n));
        end

    end

end
