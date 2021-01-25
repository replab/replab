classdef PhasedMatrixPartition < replab.Obj
% Describes a generalized partition of the cells of a matrix in disjoint subsets
%
% A standard "partition of cells of a matrix" would be a partition of the set
% $\{ (r,c) \}$ of all row-column index pairs. This standard partition (with additional constraints)
% is known as a coherent configuration:
%
% Higman, D.G. Coherent configurations. Geom Dedicata 4, 1–32 (1975)
% `<https://doi.org/10.1007/BF00147398>_`
%
% Our generalization goes in two directions:
%
% * The union of all subsets need not be the whole matrix.
% * The cells inside a subset may differ by a complex phase, and the phase needs to satisfy $\phi^n = 1$ for some integer $n$.
%
% Given a vector of coefficients ``x`` of length ``nSubsets``, we can expand it into a matrix ``X``
% such that ``X(i,j) = exp(2i*pi*phase(i,j)/phaseOrder) .* x(subsetIndex(i,j))``
% for all indices ``i,j`` such that ``subsetIndex(i,j) ~= 0``.
%
% We also have that ``[subsets{k}(1,:), subsets{k}(2,:)] = find(subsetIndex == k)``, and that
% ``[zeroSubset(1,:), zeroSubset(2,:)] = find(subsetIndex == 0)``.
%
% During construction, we assume that subsets are non-empty.

    properties (SetAccess = protected)
        nR % (integer): Number of rows in the matrix
        nC % (integer): Number of columns in the matrix
        phaseOrder % (integer): Common order of the phases
        phase % (integer(nR,nC)): Phase
        subsetIndex % (integer(nR,nC)): Subset index for each matrix cell
        zeroSubset % (integer(2,\*)): Coordinate of zero coefficients
        subsets % (cell(1,nSubsets) of integer(2,\*)): Coordinates of each subset
    end

    methods

        function self = PhasedMatrixPartition(nR, nC, phaseOrder, phase, subsetIndex, zeroSubset, subsets)
            self.nR = nR;
            self.nC = nC;
            self.phaseOrder = phaseOrder;
            self.phase = phase;
            self.subsetIndex = subsetIndex;
            self.zeroSubset = zeroSubset;
            L = cellfun(@(ss) size(ss, 2), subsets);
            assert(all(L > 0), 'Subsets (apart from the zero subset) cannot be empty');
            self.subsets = subsets;
        end

        function n = nSubsets(self)
        % Returns the number of nonzero subsets in the partition
        %
        % Returns:
        %   integer: Number of nonzero subsets
            n = length(self.subsets);
        end

        function res = normalForm(self)
            res = self.cached('normalForm', @() self.computeNormalForm);
        end

    end

    methods (Access = protected)

        function res = computeNormalForm(self)
            n = self.nSubsets;
            nR = self.nR;
            nC = self.nC;
            ind = self.subsetIndex(:);
            po = self.phaseOrder;
            changed = false;
            % ordering of subsets: the sequence given by the first occurence of each subset index must be increasing
            [sorted, sortedToUnsorted] = unique(self.subsetIndex(:)');
            sortedToUnsorted = self.subsetIndex(sortedToUnsorted);
            if sorted(1) == 0 % we do not care about 0
                sorted = sorted(2:end);
                sortedToUnsorted = sortedToUnsorted(2:end);
            end
            assert(all(sorted == 1:n));
            phase = self.phase;
            subsetIndex = self.subsetIndex;
            subsets = self.subsets;
            if any(sortedToUnsorted(:)' ~= 1:n)
                changed = true;
                unsortedToSorted = zeros(1, n);
                unsortedToSorted(sortedToUnsorted) = 1:n;
                mask = subsetIndex > 0; % only affect nonzero indices
                subsetIndex(mask) = unsortedToSorted(subsetIndex(mask));
                subsets = subsets(sortedToUnsorted);
            end
            % check that the phase order is minimal
            p = self.phase(:);
            p = unique(p(p > 0));
            if ~isempty(p)
                d = p(1);
                for i = 2:length(p)
                    d = gcd(d, p(i));
                end
                if d ~= 1
                    changed = true;
                    phase = phase/d;
                    po = po/d;
                end
            end
            % check subset by subset
            for i = 1:n
                % check that the indices of subsets are sorted by their linear index
                ss = subsets{i};
                m = size(ss, 2);
                [ss1, I] = sortrows(ss', [2 1]);
                if any(I(:)' ~= 1:m)
                    changed = true;
                    ss = ss1';
                    subsets{i} = ss;
                end
                r1 = ss(1,1);
                c1 = ss(2,1);
                % check that the phase of the first subset element is zero
                p1 = phase(r1, c1);
                if p1 ~= 0
                    changed = true;
                    ind = ss(1,:) + (ss(2,:)-1)*nR;
                    phase(ind) = mod(po + phase(ind) - p1, po);
                end
            end
            if changed
                res = replab.equi.PhasedMatrixPartition(nR, nC, po, phase, subsetIndex, self.zeroSubset, subsets);
            else
                res = self;
            end
        end

    end

    methods % Implementations

        % Obj

        function res = ne(self, rhs)
            res = ~(self == rhs);
        end

        function res = isequal(self, rhs)
            res = self == rhs;
        end

        function res = eq(self, rhs)
            if ~isa(rhs, 'replab.equi.PhasedMatrixPartition')
                res = false;
                return
            end
            lhs = self.normalForm;
            rhs = rhs.normalForm;
            res = lhs.nR == rhs.nR && lhs.nC == lhs.nR && lhs.phaseOrder == rhs.phaseOrder && ...
                  all(all(lhs.phase == rhs.phase)) && all(all(lhs.subsetIndex == rhs.subsetIndex));
        end

        function l = laws(self)
            l = replab.equi.PhasedMatrixPartitionLaws(self);
        end

    end

    methods (Static)

        function res = intersection(P1, P2)
        % Returns the partition that corresponds to the intersection of two phased matrix partitions
        %
        % "Intersection" is to be understood as the intersection of the vector spaces of matrices
        % described by ``P1`` and ``P2``.
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
                Q1 = replab.equi.PhasedMatrixPartition(nR, nC, po, P1.phase*(po/P1.phaseOrder), P1.subsetIndex, P1.zeroSubset, P1.subsets);
                Q2 = replab.equi.PhasedMatrixPartition(nR, nC, po, P2.phase*(po/P2.phaseOrder), P2.subsetIndex, P2.zeroSubset, P2.subsets);
                res = replab.equi.PhasedMatrixPartition.intersection(Q1, Q2);
                return
            end
            m = P1.phaseOrder;
            phase = P1.phase;
            subsetIndex = P1.subsetIndex;
            subsets = P1.subsets;
            zeroSubset = P1.zeroSubset;

            % We iterate through the subsets of the second PMP, and use them to join subsets of the first PMP

            % We first merge P2.zeroSubset

            % We compute the linear indices of the zero subset in the second PMP
            ind2 = P2.zeroSubset(1,:) + nR*(P2.zeroSubset(2,:)-1);
            % We compute the subset indices in the first PMP corresponding to those coordinates
            ssInds = unique(subsetIndex(ind2));
            for ssInd = ssInds(ssInds>0)
                % then all these subsets need to be set to zero
                ss = subsets{ssInd};
                subsets{ssInd} = zeros(2, 0);
                zeroSubset = horzcat(zeroSubset, ss);
                ind = ss(1,:) + nR*(ss(2,:)-1);
                subsetIndex(ind) = 0;
            end

            % Then we merge the nonzero subsets of P2
            for ssInd2 = 1:length(P2.subsets)
                % get the corresponding subset
                ss2 = P2.subsets{ssInd2};
                % compute the linear indices
                ind2 = ss2(1,:) + nR*(ss2(2,:)-1);
                % compute the subset indices in the first PMP
                ssInds = unique(subsetIndex(ind2)); % this is sorted
                isZero = false; % not yet proven to be zero
                if ssInds(1) == 0
                    % if the indices overlap the zero subset, everything is set to zero
                    isZero = true;
                    ssInds = ssInds(2:end);
                else
                    % we iterate over all subsets in the first PMP that overlap the subset of the second PMP
                    % and harmonize the phases, checking then that we do not have an inconsistency
                    % if we have an inconsistency, the subset is zero
                    for ssInd = ssInds
                        ss = subsets{ssInd};
                        ind = ss(1,:) + nR*(ss(2,:)-1);
                        common = intersect(ind, ind2);
                        delta = mod(m + P2.phase(common(1)) - phase(common(1)), m);
                        % harmonize
                        phase(common) = mod(phase(common) + P2.phase(common(1)) - phase(common(1)), m);
                        % if there is a mismatch, all the subsets are zero
                        isZero = isZero | any(phase(common) ~= P2.phase(common));
                    end
                end
                if isZero
                    % it is all zero
                    for ssInd = ssInds
                        ss = subsets{ssInd};
                        subsets{ssInd} = zeros(2, 0);
                        zeroSubset = horzcat(zeroSubset, ss);
                        ind = ss(1,:) + nR*(ss(2,:)-1);
                        subsetIndex(ind) = 0;
                        phase(ind) = 0;
                    end
                else
                    % not zero
                    remain = ssInds(1);
                    for ssInd = ssInds(2:end)
                        ss = subsets{ssInd};
                        subsets{ssInd} = zeros(2, 0);
                        subsets{remain} = horzcat(subsets{remain}, ss);
                        ind = ss(1,:) + nR*(ss(2,:)-1);
                        subsetIndex(ind) = remain;
                    end
                end
            end
            % renumber by removing empty subsets
            empty = cellfun(@isempty, subsets);
            toIndex = cumsum(~empty);
            toIndex(empty) = 0;
            fromIndex = find(~empty);
            subsets = subsets(fromIndex);
            mask = subsetIndex > 0;
            subsetIndex(mask) = toIndex(subsetIndex(mask));
            res = replab.equi.PhasedMatrixPartition(nR, nC, m, phase, subsetIndex, zeroSubset, subsets);
        end

        function res = fromPhaseAndSubsetIndexMatrices(phaseOrder, phase, subsetIndex)
        % Creates a phased matrix partition from a matrix of phases and a matrix of subset indices
            nR = size(phase, 1);
            nC = size(phase, 2);
            n = max(subsetIndex(:));
            assert(all(unique([0;subsetIndex(:)])' == 0:n));
            [I, J] = find(subsetIndex == 0);
            zeroSubset = [I(:)'; J(:)'];
            subsets = cell(1, n);
            for i = 1:n
                [I, J] = find(subsetIndex == i);
                subsets{i} = [I(:)'; J(:)'];
            end
            res = replab.equi.PhasedMatrixPartition(nR, nC, phaseOrder, phase, subsetIndex, zeroSubset, subsets);
        end

        function [res, I] = fromIndexMatrix(indexMatrix)
        % Creates a phased matrix partition from a (possibly signed) index matrix
        %
        % The index matrix contains integers. Zeros corresponds to coefficients equal to zero,
        % i.e. cells of the matrix that do not belong to a subset. The indices can be signed,
        % in which case, indices of opposite signs corresponds to coefficients of opposite value.
        %
        % If there are negative coefficients, then the result ``phaseOrder`` will be 2, otherwise it is 1.
        %
        % Note that `.PhasedMatrixPartition` requires all indices from ``1`` to ``nSubsets` to be present
        % in the partition, in other words, ``abs(indexMatrix)`` must contain all integers from ``1`` to
        % ``nSubsets``, where ``nSubsets = max(abs(indexMatrix(:)))``.
        %
        % If some indices are not present, this method will prune the corresponding subsets. The second output
        % argument ``I`` gives, for each index in ``indexMatrix``, either the index in the returned `.PhasedMatrix`
        % or ``0`` if the index does not appear in ``indexMatrix``.
        %
        % Example:
        %   >>> im = [0 1 2; -1 0 3; -2 -3 0]; % index matrix for a generic skew symmetric matrix
        %   >>> genSkewSymmetric = replab.equi.PhasedMatrixPartition.fromIndexMatrix(im);
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
            nSubsets = length(indices);
            I = cumsum(nonEmpty);
            I(~nonEmpty) = 0;
            % renumber
            subsetIndex = zeros(nR, nC);
            subsetIndex(indexMatrix~=0) = I(abs(indexMatrix(indexMatrix~=0)));
            phase = zeros(nR, nC);
            phase(indexMatrix < 0) = 1;
            phaseOrder = 1 + any(phase(:) == 1);
            res = replab.equi.PhasedMatrixPartition.fromPhaseAndSubsetIndexMatrices(phaseOrder, phase, subsetIndex);
        end

        function res = fromGeneralizedPermutations(phaseOrder, imagesR, imagesC)
        % Creates a phased matrix partition from the action of monomial representations
        %
        % It assumes that ``imagesR`` and ``imagesC`` describe generalized permutations coming from a
        % `+replab.+perm.GeneralizedSymmetricGroup` of the given ``phaseOrder``, which leave the
        % phased matrix partition we compute invariant under joint action on the rows and columns.
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
            % subsetIndex(i,j) = -1 is interpreted as "not yet visited"
            subsetIndex = -ones(nR, nC);
            subsets = cell(1, 0); % no subsets yet
            zeroSubset = zeros(2, 0);
            next = 1; % next subset index
            todo = zeros(2, 16); % stack of row/column coordinates
            Ntodo = 0; % the stack is included in column indices (1:Ntodo)
            for i = 1:nR
                for j = 1:nC
                    if subsetIndex(i, j) == -1
                        current = [i; j];
                        % for speed, we do not resize the "todo" stack, but store the index to the last element
                        todo(:, 1) = [i; j];
                        Ntodo = 1;
                        assert(phase(i, j) == 0);
                        subsetIndex(i, j) = next;
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
                                if subsetIndex(r1, c1) == -1 % not already visited
                                    % add the point to the subset
                                    phase(r1, c1) = e;
                                    subsetIndex(r1, c1) = next;
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
                            % the subset has zero coefficients, set its status
                            zeroSubset = horzcat(zeroSubset, current);
                            subsetIndex(current(1,:)+(current(2,:)-1)*nR) = 0;
                            phase(current(1,:)+(current(2,:)-1)*nR) = 0;
                            % next still points to the next free subset index
                        else
                            % add the nonzero subset
                            subsets{1,end+1} = current;
                            next = next + 1; % next free index
                        end
                    end
                end
            end
            res = replab.equi.PhasedMatrixPartition(nR, nC, m, phase, subsetIndex, zeroSubset, subsets);
        end

    end

end
