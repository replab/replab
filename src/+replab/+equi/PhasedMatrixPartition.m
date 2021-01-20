classdef PhasedMatrixPartition < replab.Obj
% Describes a generalized partition of a matrix in disjoint blocks
%
% The generalization goes in two directions:
%
% * The union of all blocks can be a subset of the whole matrix.
% * The elements inside a block may be equal up to a phase.

    properties (SetAccess = protected)
        nRows % (integer): Number of rows
        nCols % (integer): Number of columns
        nBlocks % (integer): Number of blocks
        phaseOrder % (integer): Common order of the phases
        phase % (integer(nRows,nCols)): Phase
        block % (integer(nRows,nCols)): Block index
        zeroBlock % (integer(2,\*)): Coordinate of zero coefficients
        blocks % (cell(1,nBlocks) of integer(2,\*)): Coordinates of each block
    end

    methods

        function self = PhasedMatrixPartition(nRows, nCols, phaseOrder, phase, block, zeroBlock, blocks)
            self.nRows = nRows;
            self.nCols = nCols;
            self.nBlocks = length(blocks);
            self.phaseOrder = phaseOrder;
            self.phase = phase;
            self.block = block;
            self.zeroBlock = zeroBlock;
            self.blocks = blocks;
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
            assert(P1.nRows == P2.nRows);
            assert(P1.nCols == P2.nCols);
            nR = P1.nRows;
            nC = P1.nCols;
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

        function res = fromGeneralizedSymmetricSubgroup(group, mR, mC)
        % Creates a phased matrix partition from the action of monomial representations
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group
        %   mR (`+replab.Morphism`): Morphism to a generalized symmetric group that describes a monomial representation
        %   mC (`+replab.Morphism`): Morphism to a generalized symmetric group that describes a monomial representation
        %
        % Example:
        %   >>> g1 = [1 -2 -3 -4 -5];
        %   >>> g2 = [1 4 5 2 3];
        %   >>> g3 = [1 3 2 4 -5];
        %   >>> G = replab.SignedPermutationGroup.of(g1,g2,g3);
        %   >>> mu = replab.perm.GeneralizedSymmetricGroup.morphismFromSignedPermutationGroup(G);
        %   >>> pmp = replab.equi.PhasedMatrixPartition.fromGeneralizedSymmetricSubgroup(G, mu, mu);
            if mR.target.m ~= mC.target.m
                % if the phase order is not the same, then harmonize and restart
                nR = mR.target.n;
                nC = mC.target.n;
                % common phase order
                m = lcm(mR.target.m, mC.target.m);
                % target row and column generalized symmetric groups
                targetR1 = replab.perm.GeneralizedSymmetricGroup(nR, m);
                targetC1 = replab.perm.GeneralizedSymmetricGroup(nC, m);
                % perform the conversion of the morphisms
                mR1 = mR.andThen(mR.target.type.naturalMorphism(targetR1));
                mC1 = mC.andThen(mC.target.type.naturalMorphism(targetC1));
                % call me again
                res = replab.equi.PhasedMatrixPartition(group, mR1, mC1);
                return
            end
            m = mR.target.m; % Phase order
            nR = mR.target.n; % nRows
            nC = mC.target.n; % nCols
            nG = group.nGenerators;
            % expand the images of the generators
            piR = zeros(nG, nR); % permutation part, rows
            phiR = zeros(nG, nR); % phase part, rows
            piC = zeros(nG, nC); % permutation part, cols
            phiC = zeros(nG, nC); % phase part, cols
            for i = 1:nG
                g = group.generator(i);
                % populate row and column images
                gR = mR.imageElement(g);
                piR(i, :) = gR(1, :);
                phiR(i, :) = gR(2, :);
                gC = mC.imageElement(g);
                piC(i, :) = gC(1, :);
                phiC(i, :) = gC(2, :);
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
