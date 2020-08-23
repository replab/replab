classdef SymmetricYoungIrrep < replab.Rep
% Young's orthogonal representation of a symmetric group
%
% Each irrep corresponds to an unordered partition of n.
% E.g: 4 = 2+2 and 4 = 3+1 both generate an irrep of S_4.
%
%This is a unitary representation.
%
% This is
    properties
       conjugatePartition % integer(1,:): The conjuagte partition is the partition obtained by transposing the
       % young diagram of a partition
       partition % integer(1,:):The generating partition of n
       rowFunction % integer(:,:): The i'th row is the row function for tableaux. The k'th entry is the row index k in?
       colFunction % integer(:,:): The i'th row is the column function for tableaux. Analogous
       basisHash %replab.sym.Set: Describes the row index of the row function of a basis vector.

    end

    properties(GetAccess=protected,SetAccess=protected)
         cSum
        %This represents the sum of the first n-1 elements in the partition
        % Eg: [4 2 1] => [0 4 6]
        rangeOfParts
        %This saves the array 1:(#tableax)
        underlyingRep
        % This is an underlying RepByImages Object used to quickly find the image
    end


   methods

        function self = SymmetricYoungIrrep(group, partition,form)
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
                self.isIrreducible = true;
                self.field = 'R';
                self.group = group;
                self.partition = partition;
                self.dimension = replab.sym.findPartitions.dimension(partition);
                self.conjugatePartition = replab.sym.findPartitions.conjugatePart(partition);
                self.rangeOfParts = 1:self.dimension;
                if isequal(form,'orth')
                    self.isUnitary = 1;
                end
                self.seminormalHelper();
                self.underlyingRep = self.constructRep;
        end

        function rho = image_internal(self, g)
        % Image function used by replab to calculate the images of a permutation
        %
        % Args:
        % g (permutation): Permutation whose image is calculated
        %
        % Returns:
        % rho (integer(:,:)) Image of g
            rho = self.underlyingRep.image(g);
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
            % Args:
            % k (integer): We calculate the image of the transposition
            % generators (k-1 k)
            %
            % Returns:
            % im (integer(:,:)) Image of g
                n = self.group.domainSize;
            rowFunEq = self.rowFunction(:,k) == self.rowFunction(:,k+1);
            colFunEq = self.colFunction(:,k) == self.colFunction(:,k+1);
            rowFunLess = self.rowFunction(:,k) < self.rowFunction(:,k+1);
            nInds1 = self.rangeOfParts(~rowFunEq&~colFunEq&rowFunLess);
            nInds2 = self.basisHash.find(self.rowFunction(nInds1, replab.Permutation.transposition(n, k, k+1))');
            axDistRec = 1./(self.rowFunction(nInds1,k+1)-self.rowFunction(nInds1,k) + ...
            abs(self.colFunction(nInds1,k+1)-self.colFunction(nInds1,k)));
            m1a = sparse(nInds1,nInds1,-axDistRec,self.dimension,self.dimension);
            m1b = sparse(nInds2,nInds2,axDistRec,self.dimension,self.dimension);
            if self.isUnitary
                m1c = sparse(nInds1,nInds2,sqrt(1-axDistRec.^2),self.dimension,self.dimension);
                m1d = sparse(nInds2,nInds1,sqrt(1-axDistRec.^2),self.dimension,self.dimension);
            else
                m1c = sparse(nInds1,nInds2,1-axDistRec.^2,self.dimension,self.dimension);
                m1d = sparse(nInds2,nInds1,1,self.dimension,self.dimension);
            end
            m2 = spdiags(rowFunEq,0,self.dimension,self.dimension);
            m3 = spdiags(-1*colFunEq,0,self.dimension,self.dimension);
            im = m1a+m1b+m1c+m1d+m2+m3;
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
            rep = self.group.repByImages(self.field, self.dimension, 'preimages', gens, 'images', images);
        end

   end

end
