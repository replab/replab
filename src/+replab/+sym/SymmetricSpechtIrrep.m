classdef SymmetricSpechtIrrep < replab.Obj
% An irreducible representation of a symmetric group
%
% Each irrep corresponds to an unordered partition of n.
% E.g: ``4 = 2+2`` and ``4 = 3+1`` both generate an irrep of S_4.
%
% The entries of images in these representations consist of only 0,1, and -1.
%
% This implementation is based on
% Wiltshire-Gordon, John D.; Woo, Alexander; Zajaczkowska, Magdalena (2017).
% "Specht Polytopes and Specht Matroids", https://arxiv.org/abs/1701.05277

    properties (SetAccess = protected)
        group              % (`+replab.PermutationGroup`): Symmetric group acting on `.domainSize` elements
        domainSize         % (integer): Domain size of the symmetric group
        dimension          % (integer): Representation dimension
        conjugatePartition % (integer(1,\*)): The conjuagte partition is the partition obtained by transposing the young diagram of a partition
        partition          % (integer(1,\*)): The generating partition of n
        indepRowWords      % (integer(\*,\*)): The row words; corresponding to all standard tableaux
        indepColWords      % (integer(\*,\*)): The column words; corresponding to all standard tableaux
        basis              % (integer(\*,\*)): The submatrix of the Specht matrix whose rows and columns correspond to the row and colum words of standard tableaux.
        cSum               % (integer(1,\*)): Cumulative number of boxes before the i-th line
    end

    methods (Static)

        function noRepeat = noRepeats(A)
        % Finds whether any elements in an array repeat
        %
        % Args:
        %   A (integer(1,\*)): Row vector of integers
        % Returns:
        %   noRepeats (bool): True if there are no repeats
            noRepeat = true;
            A = sort(A);
            for l = 1:length(A)-1
                if A(l) == A(l+1)
                    noRepeat = false;
                    break;
                end
            end
        end

    end

    methods

        function self = SymmetricSpechtIrrep(group, part)
        % Constructs an irreducible representation of S_n
        %
        % Args:
        %   partition (integer(1,\*)): Partition
        %   group (`+replab.Group`): Symmetric group being represented
            if sum(part) ~= group.domainSize
                error('This is not a valid Young Diagram for this domain size.')
            end
            self.group = group;
            self.domainSize = group.domainSize;
            self.dimension = replab.sym.findPartitions.dimension(part);
            self.partition = part;
            self.conjugatePartition = replab.sym.findPartitions.conjugatePart(part);
            len = numel(self.partition);
            self.cSum = [0 cumsum(self.partition(1:len-1))];
            youngLattice = replab.sym.YoungLattice(self.partition, self.domainSize);
            [self.indepRowWords, self.indepColWords, ~] = youngLattice.generateTableaux;
            % Use Young Lattice to generate words corresponding to linearly independent columns and rows
            self.basis = self.subSpecht(self.indepRowWords);
            % This is a submatrix of the Specht matrix described in the paper, chosen to be taken from linearly
            % independent columns of the matrix. We use the standard tableaux to find these and rows columns
            % and we construct only those matrix elements rather than constructing the entire Specht matrix.
            % This gives us a basis for our action.
        end

    end

    methods (Access = protected)

        function r = computeRep(self)
            preimages = self.group.generators;
            images = cellfun(@(g) self.naiveImage(g), preimages, 'uniform', 0);
            r = self.group.repByImages('R', self.dimension, 'preimages', preimages, 'images', images);
        end

    end

    methods

        function r = rep(self)
        % Returns the irreducible representation
            r = self.cached('rep', @() self.computeRep);
        end

        function specht = subSpecht(self, rowWords)
        % Finds the Specht submatrix given a list of rows and using the linearly independent columns
        %
        % Args:
        %   rowWords (integer(\*,\*)): List of row words corresponding to the rows of the submatrix
        %
        % Returns:
        %   double(\*,\*): The Specht sumbmatrix with the given rows and the known linearly independent columns
            rowInds = zeros(1, self.dimension^2);
            colInds = zeros(1, self.dimension^2);
            values = zeros(1, self.dimension^2);
            count = 0;
            % Double loop through all words corresponding to standard tableaux, i.e. a loop over all matrix elements
            for r = 1:self.dimension
                for c = 1:self.dimension
                    rowWord = rowWords(r,:);
                    colWord = self.indepColWords(c,:);
                    possiblePerm = self.cSum(rowWord) + colWord;
                    % This will be a permutation if and only if the matrix entry is nonzero
                    ent = self.entry(possiblePerm);
                    if ~isempty(ent)
                        count = count+1;
                        rowInds(count) = r;
                        colInds(count) = c;
                        values(count) = ent;
                    end
                end
            end
            specht = sparse(nonzeros(rowInds),nonzeros(colInds),nonzeros(values));
        end

        function im = naiveImage(self,perm)
            action = self.subSpecht(self.indepRowWords(:,perm));
            % A permutation acts by rearranging the rows of the Specht
            % matrix, as described in the cited paper
            im = round(self.basis/action);
            % We now can find the image of this permutation, since we know our basis
        end

        function val = entry(self, possiblePerm)
        % Finds the entry of the Specht submatrix given a row and column number
        %
        % Args:
        %   rowNum (integer): Row index
        %   colNum (integer): column index
        %
        % Returns:
        %   double: Corresponding Specht submatrix entry if its nonzero. Empty if the entry is zero.
            if replab.sym.SymmetricSpechtIrrep.noRepeats(possiblePerm)
                % This finds if the given array is a permutation
                % It is if and only if it doesn't repeat, (this is non-trivial and
                % proven in the paper)
                val = replab.Permutation.sign(possiblePerm);
                % The entry is the sign of the permutation if it is one
            else
                val = [];
                % There is no entry if it isn't a permutation
            end
         end

    end

end