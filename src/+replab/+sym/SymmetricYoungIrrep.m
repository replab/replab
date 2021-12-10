classdef SymmetricYoungIrrep < replab.Obj
% Young's orthogonal and seminormal representations of a symmetric group
%
% Each irrep corresponds to an unordered partition of n.
% E.g: ``4 = 2+2`` and ``4 = 3+1`` both generate an irrep of S_4.
%
% The "orthogonal" representation is unitary.

    properties (SetAccess = protected)
        group              % (`+replab.PermutationGroup`): Symmetric group being represented
        dimension          % (integer): Representation dimension
        form               % ('seminormal' or 'orthogonal'): Type of representation
        conjugatePartition % (integer(1,\*)): The conjuagte partition is the partition obtained by transposing the young diagram of a partition
        partition          % (integer(1,\*)): The generating partition of n
        rowFunction        % (integer(\*,\*)): The i'th row is the row function for tableaux. The k'th entry is the row index k is in
        colFunction        % (integer(\*,\*)): The i'th row is the column function for tableaux. Analogous to the row function

       % E.g the Young Tableau 1 2 5
       %                       3 4
       % has row function [1 1 2 2 1] and column function [1 2 1 2 3]
       % These are the j and j' function in Schindler Miriam

       basisHash          % (+replab.+perm.Set): Describes the row index of the row function of a Young Tableaux.
    end

    properties (Access=protected)
        cSum % Represents the sum of the first n-1 elements in the partition. Eg: [4 2 1] => [0 4 6]
        rangeOfParts % This saves the array 1:(dimension)
    end


    methods

        function self = SymmetricYoungIrrep(group, partition, form)
        % Constructs an irreducible representation of S_n
        %
        % Args:
        %   group (`replab.Group`): Group being represented
        %   partition (integer(1,\*)): Partition
        %   form ('seminormal' or 'orthogonal'): Form
            if sum(partition) ~= group.domainSize
                error('This is not a valid Young Diagram for this domain size.')
            end
            if size(partition,2) == 1
                partition = partition';
            end
            assert(size(partition, 1) == 1);
            self.group = group;
            self.form = form;
            self.dimension = replab.sym.findPartitions.dimension(partition);
            self.partition = partition;
            self.conjugatePartition = replab.sym.findPartitions.conjugatePart(partition);
            self.rangeOfParts = 1:self.dimension;
            youngLattice = replab.sym.YoungLattice(self.partition,self.group.domainSize);
            [self.rowFunction, self.colFunction, ~] = youngLattice.generateTableaux;
            %Use Young Lattice to generate row and column functions corresponding to partitions and conjugate partitions.
            self.basisHash = replab.perm.Set(self.group.domainSize);
            self.basisHash.insert(self.rowFunction');
            % Generate words corresponding to linearly independent columns
        end

        function r = rep(self)
            r = self.cached('rep', @() self.computeRep);
        end

    end

    methods (Access = protected)

        function r = computeRep(self)
            n = self.group.domainSize;
            gens = cell(1,n-1);
            for i = 1:n-1
                gens{i} = replab.Permutation.transposition(n, i, i + 1);
            end
            images = cell(1,n-1);
            for i = 1:n-1
                images{i} = self.transImage(i);
            end
            r = self.group.repByImages('R', self.dimension, 'preimages', gens, 'images', images);
        end

        function im = transImage(self,k)
        % Image function used to calculate the images of all adjacent transposition generators
        %
        % This implements the formulas as described in section III -
        % equation (9-13) of Schindler and Miriam's decomposition paper
        %
        % Args:
        %   k (integer): We calculate the image of the transposition generators (k-1 k)
        %
        % Returns:
        %   `+replab.cyclotomic` (\*,\*): Image of g
                n = self.group.domainSize;
                d = self.dimension;
                rowFunEq = self.rowFunction(:,k) == self.rowFunction(:,k+1); % Where are the row functions equal?
                colFunEq = self.colFunction(:,k) == self.colFunction(:,k+1); % Where are the column functions equal?
                % These tell us when to use Equation 10/12
                rowFunLess = self.rowFunction(:,k) < self.rowFunction(:,k+1);  % Where are the j_k less than j_(k+1)?
                neitherAndIsR = self.rangeOfParts(~rowFunEq & ~colFunEq & rowFunLess);
                transp = replab.Permutation.transposition(n, k, k+1);
                neitherAndIsRPrime = self.basisHash.find(self.rowFunction(neitherAndIsR, transp)');
                % These tell us when to use Equation 11/13 and which are r and r', as described in the paper.
                axDistRec = (self.rowFunction(neitherAndIsR,k+1)-self.rowFunction(neitherAndIsR,k) + ...
                    abs(self.colFunction(neitherAndIsR,k+1)-self.colFunction(neitherAndIsR,k)));
                axDistRec = replab.cyclotomic.eye(1) ./ axDistRec;
                % This is the axial distance described in equation 9 but
                % only calculated where we need it (for equation 11/13)
                m1a = replab.cyclotomic.sparse(neitherAndIsR, neitherAndIsR, -axDistRec, d, d);
                m1b = replab.cyclotomic.sparse(neitherAndIsRPrime, neitherAndIsRPrime, axDistRec, d, d);
                % The m1's describes the matrix elements from Equation 11/13.
                % m1a and m1b describe the common matrix elements in Eq 11 and
                % 13
                axDistRec2 = axDistRec.*axDistRec;
                if strcmp(self.form, 'orthogonal') % preset this if we are in the orthogonal representation
                    m1c = replab.cyclotomic.sparse(neitherAndIsR, neitherAndIsRPrime, sqrt(1-axDistRec2), d, d);
                    m1d = replab.cyclotomic.sparse(neitherAndIsRPrime, neitherAndIsR, sqrt(1-axDistRec2), d, d);
                    % These are the unqiue matrix elements from Eq 13
                else
                    m1c = replab.cyclotomic.sparse(neitherAndIsR, neitherAndIsRPrime, 1-axDistRec2, d, d);
                    m1d = replab.cyclotomic.sparse(neitherAndIsRPrime, neitherAndIsR, ones(1, length(neitherAndIsR)), d, d);
                    % These are the unique matrix elements from Eq 11
                end
                m2 = double(full(spdiags(rowFunEq, 0, d, d)));
                m3 = full(spdiags(-1*colFunEq, 0, d, d));
                % These are the matrix elements from Eq 10 and 12 (note
                % that they are the same)
                im = m1a+m1b+m1c+m1d+m2+m3;
                % Add all of the matrix elements together
        end

   end

end
