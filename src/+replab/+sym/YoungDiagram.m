classdef YoungDiagram < replab.Str

    properties (SetAccess = protected)
        n % integer: Number of boxes in the Young diagram
        partition % (integer(1,\*)): Partition corresponding to the Young diagram
    end

    methods (Static)

        function I = allPartitions(n, m)
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
                I = horzcat(I, cellfun(@(part) [i part], replab.sym.YoungDiagram.allPartitions(n - i, i), 'uniform', 0));
            end
        end

        function Y = allYoungDiagrams(n)
        % Returns all Young Diagrams for the symmetric group of given degree
        %
        % Args:
        %   n (integer): Symmetric group domain size
        %
        % Returns:
        %   cell(1,\*) of `.YoungDiagram`: Young diagrams
            Y = cellfun(@(p) replab.sym.YoungDiagram(p), replab.sym.YoungDiagram.allPartitions(n), 'uniform', 0);
        end

        function [gens order] = centralizerBlocks(n, start, blockSize, nBlocks)
            order = vpi(blockSize)^nBlocks*factorial(vpi(nBlocks));
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

    methods

        function self = YoungDiagram(partition)
            self.partition = partition;
            self.n = sum(partition);
        end

        function C = conjugacyClass(self)
        % Returns the conjugacy class of the symmetric group corresponding to the Young Diagram
            rep = self.minLexPermutation;
            repCent = self.centralizerMinLexPermutation;
            C = replab.ConjugacyClass(replab.S(self.n), rep, repCent);
        end

        function p = minLexPermutation(self)
        % Returns the minimal permutation under lexicographic ordering that has its cycle structure given by the partition of this Young diagram
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
            order = vpi(1);
            while j <= length(part) % for each part size
                j1 = find([part(j+1:end) 0] ~= part(j), 1) + j; % find when the next part size starts
                blockSize = part(j);
                nBlocks = j1 - j;
                % computes the wreath product of C(blockSize) by S(nBlocks)
                [gensJ orderJ] = replab.sym.YoungDiagram.centralizerBlocks(self.n, ind, blockSize, nBlocks);
                % append generators and multiply order
                gens = horzcat(gens, gensJ);
                order = order * orderJ;
                ind = ind + nBlocks * blockSize;
                j = j1;
            end
            % construct the group
            G = replab.PermutationGroup(self.n, gens, order, replab.SymmetricGroup.make(self.n));
        end


        function [above, left] = aboveLeft(self)
        % Calculates positional information of the Young diagram
        %
        % Returns
        % -------
        % above:
        %   integer(1,\*): Index of the box immediately to the top in the Young diagram, 0 if none
        % left:
        %   integer(1,\*): Index of the box immediately on the left in the Young diagram 0 if none
            part = self.partition;
            n = sum(part);
            m = numel(part);
            cSum = [0 cumsum(part(1:m-1))];
            above = zeros(1, n);
            left = 0:(n-1);
            left(cSum + 1)=0;
            for j = 2:m
                inds = 1:part(j);
                above(cSum(j) + inds) = cSum(j-1) + inds;
            end
        end

    end

end
