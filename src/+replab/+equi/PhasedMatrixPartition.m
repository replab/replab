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

        function self = PhasedMatrixPartition(nRows, nCols, nBlocks, phaseOrder, phase, block, zeroBlock, blocks)
            self.nRows = nRows;
            self.nCols = nCols;
            self.nBlocks = nBlocks;
            self.phaseOrder = phaseOrder;
            self.phase = phase;
            self.block = block;
            self.zeroBlock = zeroBlock;
            self.blocks = blocks;
        end

    end

    methods (Static)

        function res = fromGeneralizedSymmetricSubgroup(group, mR, mC)
        % Creates a phased matrix partition from the action of monomial representations
            if mR.target.m ~= mC.target.m
                % if the phase order is not the same, then harmonize and restart
                nR = mR.target.n;
                nC = mC.target.n;
                m = lcm(mR.target.m, mC.target.m);
                targetR1 = replab.perm.GeneralizedSymmetricGroup(nR, m);
                targetC1 = replab.perm.GeneralizedSymmetricGroup(nC, m);
                mR1 = mR.andThen(mR.target.type.naturalMorphism(targetR1));
                mC1 = mC.andThen(mC.target.type.naturalMorphism(targetC1));
                res = replab.equi.PhasedMatrixPartition(group, mR1, mC1);
            else
                m = mR.target.m; % Phase order
                nR = mR.target.n; % nRows
                nC = mC.target.n; % nCols
                nG = group.nGenerators;
                piR = zeros(nG, nR);
                phiR = zeros(nG, nR);
                piC = zeros(nG, nC);
                phiC = zeros(nG, nC);
                for i = 1:nG
                    g = group.generator(i);
                    gR = mR.imageElement(g);
                    piR(i, :) = gR(1, :);
                    phiR(i, :) = gR(2, :);
                    gC = mC.imageElement(g);
                    piC(i, :) = gC(1, :);
                    phiC(i, :) = gC(2, :);
                end
                phase = zeros(nR, nC);
                block = -ones(nR, nC);
                blocks = cell(1, 0);
                zeroBlock = zeros(2, 0);
                next = 1;
                todo = zeros(2, 16);
                Ntodo = 0;
                for i = 1:nR
                    for j = 1:nC
                        if block(i, j) == -1
                            current = [i; j];
                            todo(:, 1) = [i; j];
                            Ntodo = 1;
                            assert(phase(i, j) == 0);
                            block(i, j) = next;
                            isZero = false;
                            while Ntodo > 0
                                r = todo(1, Ntodo);
                                c = todo(2, Ntodo);
                                Ntodo = Ntodo - 1;
                                for k = 1:nG
                                    r1 = piR(k, r);
                                    c1 = piC(k, c);
                                    e = mod(phase(r, c) + phiR(k, r) - phiC(k, c) + m, m);
                                    if block(r1, c1) == -1
                                        phase(r1, c1) = e;
                                        block(r1, c1) = next;
                                        Ntodo = Ntodo + 1;
                                        todo(:, Ntodo) = [r1; c1];
                                        current(:, end+1) = [r1; c1];
                                    elseif ~isZero
                                        if phase(r1, c1) ~= e
                                            isZero = true;
                                        end
                                    end
                                end
                            end
                            if isZero
                                zeroBlock = horzcat(zeroBlock, current);
                                block(current(1,:)+(current(2,:)-1)*nR) = 0;
                                phase(current(1,:)+(current(2,:)-1)*nR) = 0;
                            else
                                blocks{1,end+1} = current;
                                next = next + 1;
                            end
                        end
                    end
                end
                res = replab.equi.PhasedMatrixPartition(nR, nC, next - 1, m, phase, block, zeroBlock, blocks);

            end
        end % fromGeneralizedSymmetricSubgroup

    end

end
