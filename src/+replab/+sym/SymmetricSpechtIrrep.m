classdef SymmetricSpechtIrrep < replab.Rep

% An irreducible representation of a symmetric group
%
% Each irrep corresponds to an unordered partition of n.
% E.g: 4 = 2+2 and 4 = 3+1 both generate an irrep of S_4.
%
%The entries of images in these represntations consist of only 0,1, and -1.
%
% This implementation is based on
% Wiltshire-Gordon, John D.; Woo, Alexander; Zajaczkowska, Magdalena (2017).
% "Specht Polytopes and Specht Matroids"
% https://arxiv.org/abs/1701.05277

    properties
        conjugatePartition % integer(1,:): The conjuagte partition is the partition obtained by transposing the
                           % young diagram of a partition
        partition % integer(1,:):The generating partition of n
        indepRowWords % integer(:,:): The row words; corresponding to all standard tableaux
        indepColWords% integer(:,:): The column words; corresponding to all standard tableaux
        basis % integer(:,:): The submatrix of the Specht matrix whose rows and columns correspond to the row and colum
              %words of standard tableaux.
    end

    properties(GetAccess=protected,SetAccess=protected)
        cSum
        underlyingRep
        % This is an underlying RepByImages object used to quickly find the image
    end

    methods

        function self = SymmetricSpechtIrrep(group, part)
        % Constructs an irreducible representation of S_n
        %
        % Args:
        % partition (integer(1,:)): Partition
        % group (`replab.Group`): Symmetric group being represented
        %
        % Returns:
        % self (`replab.Rep`): Specht representation of the group
            if sum(part) ~= group.domainSize
                error('This is not a valid Young Diagram for this domain size.')
            end
            self.field = 'R';
            self.group = group;
            self.partition = part;
            self.dimension = replab.sym.findPartitions.dimension(part);
            self.conjugatePartition = replab.sym.findPartitions.conjugatePart(part);
            symIrrepHelper(self);
            self.underlyingRep = replab.RepByImages.fromImageFunction(self.group, self.field, self.dimension, @(g) naiveImage(self,g));

        end

        function rho = image_internal(self, g)
        % Image function used by replab to calculate the images of a permutation
        %
        % Args:
        % g (permutation): Permutation whose image is calculated
        %
        % Returns:
        % rho (integer(:,:)) Image of g
            rho = round(self.underlyingRep.image(g));
        end

    end

    methods(Access = protected)

        function symIrrepHelper(self)
        % Helper function for constructor
            len = numel(self.partition);
            self.cSum = [0 cumsum(self.partition(1:len-1))];
            youngLattice = replab.sym.YoungLattice(self.partition,self.group.domainSize);
            [self.indepRowWords,self.indepColWords,~] = youngLattice.generateTableaux;
            %Use Young Lattice to generate words corresponding to linearly independent
            %columns and rows
            self.basis = self.subSpecht(self.indepRowWords);
        end

        function specht = subSpecht(self,rowWords)
        % Finds the Specht submatrix given a list of rows
        % and using the linearly independent columns
        %
        % Args:
        % rowWords (integer(:.:)): List of row words
        % corresponding to the rows of the submatrix
        %
        % Returns:
        % specht (double(:,:)): The specht sumbmatrix with
        % the given rows and the known linearly independent
        % columns
            rowInds = zeros(1,self.dimension^2);
            colInds = zeros(1,self.dimension^2);
            values = zeros(1,self.dimension^2);
            count = 0;
            for r = 1:self.dimension
                for c = 1:self.dimension
                    ent = entry(r,c);
                    if ~isempty(ent)
                        count = count+1;
                        rowInds(count) = r;
                        colInds(count) = c;
                        values(count) = ent;
                    end
                end
            end
            specht = sparse(nonzeros(rowInds),nonzeros(colInds),nonzeros(values));

            function val = entry(rowNum,colNum)
            % Finds the entry of the Specht submatrix given a row
            % and column number
            %
            % Args:
            % rowNum (integer): Row index
            % colNum (integer): column index
            %
            % Returns:
            % val (double): Corresponding Specht submatrix entry if
            % its nonzero. Empty if the entry is zero.
                row = rowWords(rowNum,:);
                col = self.indepColWords(colNum,:);
                perm = self.cSum(row) + col;
                if noRepeats(perm)
                    val = PermSign(perm);
                else
                    val = [];
                end


                function noRepeat = noRepeats(A)
                % Finds whether any elements in an array repeat
                %
                % Args:
                % A (integer(1,:)):
                % Returns:
                % noRepeats (bool): True if there are no repeats
                    noRepeat = true;
                    A = sort(A);
                    for l=1:(self.group.domainSize-1)
                        if A(l)==A(l+1)
                            noRepeat = false;
                            break;
                        end
                    end
                end

                function sign = PermSign(perm)
                %Returns sign of a given permuatation
                %perm (row vector): Vector representing a permutation (e.g. [3 2 1 4])
                %sign (float): sign of the permutation
                    x = perm;
                    n = length(x);
                    oddOrEven = 0; %Records whether the permutation is odd or even
                    for k = 1:n
                        if x(k) == 0 || x(k) == k %Skip over one cycles and numbers that have been cycled through
                            continue
                        end
                        cycleSize = -1; %The first element in a cycle isn't counted
                        q = k;
                        while x(q) ~= 0
                            pHold = x(q);
                            x(q) = 0;
                            q = pHold;
                            cycleSize = cycleSize + 1;
                        end
                        if cycleSize > 0
                            oddOrEven = oddOrEven + cycleSize; %At the end, this will match the parity (even/odd) of the permuation
                        end
                    end
                    sign = (-1)^mod(round(oddOrEven),2); %Sign of permutation
                end
            end
        end

        function im = naiveImage(self,perm)
            action = self.subSpecht(self.indepRowWords(:,perm));
            im = round(self.basis/action);
        end

    end

end