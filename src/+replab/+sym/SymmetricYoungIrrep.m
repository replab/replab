classdef SymmetricYoungIrrep < replab.Rep
% Young's orthogonal and seminormal representations of a symmetric group
%
% Each irrep corresponds to an unordered partition of n.
% E.g: 4 = 2+2 and 4 = 3+1 both generate an irrep of S_4.
%
%This is a unitary representation.
%
% This is
    properties
       conjugatePartition % integer(1,:): The conjuagte partition is the partition obtained by transposing the young diagram of a partition
       partition % integer(1,:):The generating partition of n
       rowFunction % integer(:,:): The i'th row is the row function for tableaux. The k'th entry is the row index k is in
       colFunction % integer(:,:): The i'th row is the column function for tableaux. Analogous to the row function

       % E.g the Young Tableax 1 2 5
       %                       3 4
       % has row function [1 1 2 2 1] and column function [1 2 1 2 3]
       % These are the j and j' function in Schindler Miriam

       basisHash % replab.sym.Set: Describes the row index of the row function of a Young Tableaux.
    end

    properties(GetAccess=protected,SetAccess=protected)
        cSum % Represents the sum of the first n-1 elements in the partition. Eg: [4 2 1] => [0 4 6]
        rangeOfParts % This saves the array 1:(dimension)
        underlyingRep % This is an underlying RepByImages Object used to quickly find the image
    end


    methods

        function self = SymmetricYoungIrrep(group, partition, form)
        % Constructs an irreducible representation of S_n
        %
        % Args:
        % partition (integer(1,:)): Partition
        % group (`replab.Group`): Group being representation
            if sum(partition) ~= group.domainSize
                error('This is not a valid Young Diagram for this domain size.')
            end
            if size(partition,2) == 1
                partition = partition';
            end
            assert(size(partition,1) == 1);
            dimension = replab.sym.findPartitions.dimension(partition);
            self@replab.Rep(group, 'R', dimension, 'isIrreducible', true, 'isUnitary', strcmp(form, 'orth'));
            self.partition = partition;
            self.conjugatePartition = replab.sym.findPartitions.conjugatePart(partition);
            self.rangeOfParts = 1:self.dimension;
            self.seminormalHelper();
            self.underlyingRep = self.constructRep;
        end

    end

    methods

        function b = canComputeType(self, type)
            b = self.underlyingRep.canComputeType(type);
        end

        function rho = image(self, g, varargin)
        % Image function used by replab to calculate the images of a permutation
        %
        % Args:
        % g (permutation): Permutation whose image is calculated
        %
        % Returns:
        % rho (integer(:,:)) Image of g
            rho = self.underlyingRep.image(g, varargin{:});
        end

    end

    methods (Access = protected)

        function e = computeErrorBound(self)
            e = self.underlyingRep.errorBound;
        end

    end

    methods (Access = protected)

        function seminormalHelper(self)
        % Helper function for constructor

            youngLattice = replab.sym.YoungLattice(self.partition,self.group.domainSize);
            [self.rowFunction,self.colFunction,~] = youngLattice.generateTableaux;
            %Use Young Lattice to generate row and
            %column functions corresponding to partitions and
            %conjugate partitions.
            self.basisHash = replab.perm.Set(self.group.domainSize);
            self.basisHash.insert(self.rowFunction');
            %Generate words corresponding to linearly independent columns
        end


        function im = transImage(self,k)
            % Image function used to calculate the images of all adjacent transposition generators
            %
            % This implements the formulas as described in section III -
            % equation (9-13) of Schindler and Miriam's decomposition paper
            %
            % Args:
            % k (integer): We calculate the image of the transposition
            % generators (k-1 k)
            %
            % Returns:
            % im (integer(:,:)) Image of g
                n = self.group.domainSize;
                rowFunEq = self.rowFunction(:,k) == self.rowFunction(:,k+1); %Where are the row functions equal?
                colFunEq = self.colFunction(:,k) == self.colFunction(:,k+1); %Where are the column functions equal?
                %These tell us when to use Equation 10/12
                rowFunLess = self.rowFunction(:,k) < self.rowFunction(:,k+1);  %Where are the j_k less than j_(k+1)?
                neitherAndIsR = self.rangeOfParts(~rowFunEq&~colFunEq&rowFunLess);
                neitherAndIsRPrime = self.basisHash.find(self.rowFunction(neitherAndIsR, ...
                    replab.Permutation.transposition(n, k, k+1))');
                %These tell us when to use Equation 11/13 and which are r
                %and r', as described in the paper.
                axDistRec = 1./(self.rowFunction(neitherAndIsR,k+1)-self.rowFunction(neitherAndIsR,k) + ...
                    abs(self.colFunction(neitherAndIsR,k+1)-self.colFunction(neitherAndIsR,k)));
                % This is the axial distance described in equation 9 but
                % only calculated where we need it (for equation 11/13)
                m1a = sparse(neitherAndIsR,neitherAndIsR,-axDistRec,self.dimension,self.dimension);
                m1b = sparse(neitherAndIsRPrime,neitherAndIsRPrime,axDistRec,self.dimension,self.dimension);
                % The m1's describes the matrix elements from Equation 11/13.
                % m1a and m1b describe the common matrix elements in Eq 11 and
                % 13
                if self.isUnitary %preset this if we are in the orthogonal representation
                    m1c = sparse(neitherAndIsR,neitherAndIsRPrime,sqrt(1-axDistRec.^2),self.dimension,self.dimension);
                    m1d = sparse(neitherAndIsRPrime,neitherAndIsR,sqrt(1-axDistRec.^2),self.dimension,self.dimension);
                    % These are the unqiue matrix elements from Eq 13
                else
                    m1c = sparse(neitherAndIsR,neitherAndIsRPrime,1-axDistRec.^2,self.dimension,self.dimension);
                    m1d = sparse(neitherAndIsRPrime,neitherAndIsR,1,self.dimension,self.dimension);
                    % These are the unique matrix elements from Eq 11
                end
                m2 = spdiags(rowFunEq,0,self.dimension,self.dimension);
                m3 = spdiags(-1*colFunEq,0,self.dimension,self.dimension);
                % These are the matrix elements from Eq 10 and 12 (note
                % that they are the same)
                im = m1a+m1b+m1c+m1d+m2+m3;
                % Add all of the matrix elements together
        end

        function rep = constructRep(self)
            n = self.group.domainSize;
            gens = cell(1,n-1);
            for i = 1:n-1
                gens{i} = replab.Permutation.transposition(n, i, i + 1);
            end
            images = cell(1,n-1);
            for i = 1:n-1
                images{i} = self.transImage(i);
            end
            rep = self.group.repByImages(self.field, self.dimension, 'preimages', gens, 'images', images, 'isUnitary', self.isUnitary);
        end

   end

end
