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
            else
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
