classdef PhasedMatrixPartition < replab.Obj
% Describes a generalized partition of a matrix in disjoint blocks
%
% The generalization goes in two directions:
%
% * The union of all blocks need not be the whole matrix.
% * The elements inside a block may differ by a complex phase, which needs to be a rational root of unity.
%
% For a vector of coefficients ``x`` of length ``nBlocks``, it describes a matrix ``X`` such that
% ``X(i,j) = exp(2i*pi*phase(i,j)/phaseOrder) .* x(block(i,j))`` for the indices ``i,j`` such that
% ``block(i,j) > 0``.
%
% We also have that ``[blocks{k}(1,:), blocks{k}(2,:)] = find(block == k)``, and that
% ``[zeroBlock(1,:), zeroBlock(2,:)] = find(block == 0)``.
%
% During construction, we assume that blocks are non-empty.

    properties (SetAccess = protected)
        nR % (integer): Number of rows in the matrix
        nC % (integer): Number of columns in the matrix
        phaseOrder % (integer): Common order of the phases
        phase % (integer(nR,nC)): Phase
        block % (integer(nR,nC)): Block index
        zeroBlock % (integer(2,\*)): Coordinate of zero coefficients
        blocks % (cell(1,nBlocks) of integer(2,\*)): Coordinates of each block
    end

    methods

        function self = PhasedMatrixPartition(nR, nC, phaseOrder, phase, block, zeroBlock, blocks)
            self.nR = nR;
            self.nC = nC;
            self.phaseOrder = phaseOrder;
            self.phase = phase;
            self.block = block;
            self.zeroBlock = zeroBlock;
            L = cellfun(@(blk) size(blk, 2), blocks);
            assert(all(L > 0), 'Blocks cannot be empty');
            self.blocks = blocks;
        end

        function n = nBlocks(self)
        % Returns the number of nonzero blocks in the partition
        %
        % Returns:
        %   integer: Number of nonzero blocks
            n = length(self.blocks);
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.equi.PhasedMatrixPartitionLaws(self);
        end

    end

    methods (Static)

        function res = intersection(P1, P2)
        % Returns the partition that corresponds to the intersection of two phased matrix partitions
        %
        % Args:
        %   P1 (`.PhasedMatrixPartition`): First partition
        %   P2 (`.PhasedMatrixPartition`): Second partition
        %
        % Returns:
        %   `.PhasedMatrixPartition`: Intersection
            assert(P1.nR == P2.nR);
            assert(P1.nC == P2.nC);
            nR = P1.nR;
            nC = P1.nC;
            if P1.phaseOrder ~= P2.phaseOrder
                % if the phase order is not the same, then harmonize and restart
                po = lcm(P1.phaseOrder, P2.phaseOrder);
                Q1 = replab.equi.PhasedMatrixPartition(nR, nC, P1.nBlocks, po, P1.phase*(po/P1.phaseOrder), P1.block, P1.zeroBlock, P1.blocks);
                Q2 = replab.equi.PhasedMatrixPartition(nR, nC, P2.nBlocks, po, P2.phase*(po/P2.phaseOrder), P2.block, P2.zeroBlock, P2.blocks);
                res = replab.equi.PhasedMatrixPartition.intersection(Q1, Q2);
                return
            end
            m = P1.phaseOrder;
            phase = P1.phase;
            block = P1.block;
            blocks = P1.blocks;
            zeroBlock = P1.zeroBlock;

            % We iterate through the blocks of the second PMP, and use them to join blocks of the first PMP

            % We first merge P2.zeroBlock

            % We compute the linear indices of the zero block in the second PMP
            ind2 = P2.zeroBlock(1,:) + nR*(P2.zeroBlock(2,:)-1);
            % We compute the block indices in the first PMP corresponding to those coordinates
            blkInds = unique(block(ind2));
            for blkInd = blkInds(blkInds>0)
                % then all these blocks need to be set to zero
                blk = blocks{blkInd};
                blocks{blkInd} = zeros(2, 0);
                zeroBlock = horzcat(zeroBlock, blk);
                ind = blk(1,:) + nR*(blk(2,:)-1);
                block(ind) = 0;
            end

            % Then we merge the nonzero blocks of P2
            for blkInd2 = 1:length(P2.blocks)
                % get the corresponding block
                blk2 = P2.blocks{blkInd2};
                % compute the linear indices
                ind2 = blk2(1,:) + nR*(blk2(2,:)-1);
                % compute the block indices in the first PMP
                blkInds = unique(block(ind2)); % this is sorted
                isZero = false; % not yet proven to be zero
                if blkInds(1) == 0
                    % if the indices overlap the zero block, everything is set to zero
                    isZero = true;
                    blkInds = blkInds(2:end);
                else
                    % we iterate over all blocks in the first PMP that overlap the block of the second PMP
                    % and harmonize the phases, checking then that we do not have an inconsistency
                    % if we have an inconsistency, the block is zero
                    for blkInd = blkInds
                        blk = blocks{blkInd};
                        ind = blk(1,:) + nR*(blk(2,:)-1);
                        common = intersect(ind, ind2);
                        delta = mod(m + P2.phase(common(1)) - phase(common(1)), m);
                        % harmonize
                        phase(common) = mod(phase(common) + P2.phase(common(1)) - phase(common(1)), m);
                        % if there is a mismatch, all the blocks are zero
                        isZero = isZero | any(phase(common) ~= P2.phase(common));
                    end
                end
                if isZero
                    % it is all zero
                    for blkInd1 = blkInds1
                        blk = blocks{blkInd1};
                        blocks{blkInd1} = zeros(2, 0);
                        zeroBlock = horzcat(zeroBlock, blk);
                        ind1 = blk(1,:) + nR*(blk(2,:)-1);
                        block(ind1) = 0;
                        phase(ind1) = 0;
                    end
                else
                    % not zero
                    remain = blkInds(1);
                    for blkInd1 = blkInds(2:end)
                        blk = blocks{blkInd1};
                        blocks{blkInd1} = zeros(2, 0);
                        blocks{remain} = horzcat(blocks{remain}, blk);
                        ind1 = blk(1,:) + nR*(blk(2,:)-1);
                        block(ind1) = remain;
                    end
                end
            end
            % renumber by removing empty blocks
            empty = cellfun(@isempty, blocks);
            toIndex = cumsum(~empty);
            toIndex(empty) = 0;
            fromIndex = find(~empty);
            blocks = blocks(fromIndex);
            block(block>0) = toIndex(block(block>0));
            res = replab.equi.PhasedMatrixPartition(nR, nC, m, phase, block, zeroBlock, blocks);
        end

        function [res, I] = fromIndexMatrix(indexMatrix)
        % Creates a phased matrix partition from a (possibly signed) index matrix
        %
        % The index matrix contains integers. Zeros corresponds to coefficients equal to zero,
        % i.e. cells of the matrix that do not belong to a block. The indices can be signed,
        % in which case, indices of opposite signs corresponds to coefficients of opposite value
        %
        % If there are negative coefficients, then the result ``phaseOrder`` will be 2, otherwise it is 1.
        %
        % Note that `.PhasedMatrixPartition` requires all indices from ``1`` to ``nBlocks` to be present
        % in the partition, in other words, ``abs(indexMatrix)`` must contain all integers from ``1`` to
        % ``nBlocks``, where ``nBlocks = max(abs(indexMatrix(:)))``.
        %
        % If some indices are not present, this method will prune the corresponding blocks. The second output
        % argument ``I`` gives, for each index in ``indexMatrix``, either the index in the returned `.PhasedMatrix`
        % or ``0`` if the index does not appear in ``indexMatrix``.
        %
        % Args:
        %   indexMatrix (integer(\*,\*)): Matrix of coefficient indices, possibly signed
        %
        % Returns
        % -------
        %   res: `.PhasedMatrixPartition`
        %     The corresponding matrix partition
        %   I: integer(1,\*)
        %     Index renumbering scheme
            nR = size(indexMatrix, 1);
            nC = size(indexMatrix, 2);
            % compute index reordering
            indices = unique(abs(indexMatrix(:)));
            indices = indices(indices > 0);
            n = max(indices);
            nonEmpty = false(1, n);
            nonEmpty(indices) = true;
            nBlocks = length(indices);
            I = cumsum(nonEmpty);
            I(~nonEmpty) = 0;
            % renumber
            block = zeros(nR, nC);
            block(indexMatrix~=0) = I(abs(indexMatrix(indexMatrix~=0)));
            phase = zeros(nR, nC);
            phase(indexMatrix < 0) = 1;
            phaseOrder = 1 + any(phase(:) == 1);
            % compute blocks
            [I, J] = find(block == 0);
            zeroBlock = [I(:)'; J(:)'];
            blocks = cell(1, nBlocks);
            for i = 1:nBlocks
                [I, J] = find(block == i);
                blocks{i} = [I(:)'; J(:)'];
            end
            res = replab.equi.PhasedMatrixPartition(nR, nC, phaseOrder, phase, block, zeroBlock, blocks);
        end

        function res = fromGeneralizedPermutations(phaseOrder, imagesR, imagesC)
        % Creates a phased matrix partition from the action of monomial representations
        %
        % Args:
        %   phaseOrder (integer): Phase order
        %   imagesR (cell(1,\*) of generalized permutations): Generalized permutations acting on the row space (non-empty)
        %   imagesC (cell(1,\*) of generalized permutations): Generalized permutatioms acting on the column space (non-empty)
        %
        % Example:
        %   >>> g1 = [1 -2 -3 -4 -5];
        %   >>> g2 = [1 4 5 2 3];
        %   >>> g3 = [1 3 2 4 -5];
        %   >>> G = replab.SignedPermutationGroup.of(g1,g2,g3);
        %   >>> mu = replab.perm.GeneralizedSymmetricGroup.morphismFromSignedPermutationGroup(G);
        %   >>> h1 = mu.imageElement(g1); h2 = mu.imageElement(g2); h3 = mu.imageElement(g3);
        %   >>> pmp = replab.equi.PhasedMatrixPartition.fromGeneralizedPermutations(2, {h1 h2 h3}, {h1 h2 h3});
            m = phaseOrder; % Phase order
            nR = size(imagesR{1}, 2);
            nC = size(imagesC{1}, 2);
            nG = length(imagesR);
            assert(nG == length(imagesC));
            % expand the images of the generators
            piR = zeros(nG, nR); % permutation part, rows
            phiR = zeros(nG, nR); % phase part, rows
            piC = zeros(nG, nC); % permutation part, cols
            phiC = zeros(nG, nC); % phase part, cols
            for i = 1:nG
                iR = imagesR{i};
                iC = imagesC{i};
                piR(i,:) = iR(1,:);
                phiR(i,:) = iR(2,:);
                piC(i,:) = iC(1,:);
                phiC(i,:) = iC(2,:);
            end
            phase = zeros(nR, nC);
            % block(i,j) = -1 is interpreted as "not yet visited"
            block = -ones(nR, nC);
            blocks = cell(1, 0); % no blocks yet
            zeroBlock = zeros(2, 0);
            next = 1; % next block index
            todo = zeros(2, 16); % stack of row/column coordinates
            Ntodo = 0; % the stack is included in column indices (1:Ntodo)
            for i = 1:nR
                for j = 1:nC
                    if block(i, j) == -1
                        current = [i; j];
                        % for speed, we do not resize the "todo" stack, but store the index to the last element
                        todo(:, 1) = [i; j];
                        Ntodo = 1;
                        assert(phase(i, j) == 0);
                        block(i, j) = next;
                        isZero = false;
                        % visit all points whose images need to be checked in turn
                        while Ntodo > 0
                            % pop the last element from the stack
                            r = todo(1, Ntodo);
                            c = todo(2, Ntodo);
                            Ntodo = Ntodo - 1;
                            % compute the image under all generators
                            for k = 1:nG
                                r1 = piR(k, r);
                                c1 = piC(k, c);
                                e = mod(phase(r, c) + phiR(k, r) - phiC(k, c) + m, m);
                                if block(r1, c1) == -1 % not already visited
                                    % add the point to the block
                                    phase(r1, c1) = e;
                                    block(r1, c1) = next;
                                    current(:, end+1) = [r1; c1];
                                    % add the point to the list of points to check
                                    Ntodo = Ntodo + 1;
                                    todo(:, Ntodo) = [r1; c1];
                                elseif ~isZero && phase(r1, c1) ~= e
                                    isZero = true;
                                end
                            end
                        end
                        if isZero
                            % the block has zero coefficients, set its status
                            zeroBlock = horzcat(zeroBlock, current);
                            block(current(1,:)+(current(2,:)-1)*nR) = 0;
                            phase(current(1,:)+(current(2,:)-1)*nR) = 0;
                            % next still points to the next free block index
                        else
                            % add the nonzero block
                            blocks{1,end+1} = current;
                            next = next + 1; % next free index
                        end
                    end
                end
            end
            res = replab.equi.PhasedMatrixPartition(nR, nC, m, phase, block, zeroBlock, blocks);
        end

    end

end
