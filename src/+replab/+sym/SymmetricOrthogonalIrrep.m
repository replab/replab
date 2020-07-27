classdef SymmetricOrthogonalIrrep < replab.Rep
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
       basisOrder %struct: Describes the row a partition is in.

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
        function self = SymmetricOrthogonalIrrep(group, partition)
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
                self.dimension = replab.sym.partition.dimension(partition);
                self.conjugatePartition = replab.sym.partition.conjugatePart(partition);
                self.rangeOfParts = 1:self.dimension;
                self.isUnitary = 1;
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
   methods(Access =protected)
       
        function seminormalHelper(self)
        % Helper function for constructor
            baseWords = replab.sym.words(self.partition,self.conjugatePartition,'char'); 
            %Generate words corresponding to partition and
            %conjuagate partition.
            self.basisOrder = struct;
            [self.rowFunction,self.colFunction] = tableaux;  
            %Generate words corresponding to linearly independent columns
            function [jFun,jPrimeFun] = tableaux
            % Enumerate standard Young tableaux for a given partition
            %
            % Indices in the Young diagram increase left to right then top to bottom
            %
            % These tableaux are in one-to-one correspondence with rearrangements of words
            %
            % Returns:
            %  rowWords (integer(:,:): A matrix whose rows enumerate all row words corresponding to standard Young tableaux
            %  colWords (integer(:,:): A matrix whose rows enumerate all column words corresponding to standard Young tableaux
                n = sum(self.partition);
                %stores all indices and values for sparse matrix construction,
                %The ith row stores values for the the (i)-(i+1) transpoition
                %image
                rowCount = 0;
                [left_index, top_index] = aboveLeft(self.partition);
                rec(zeros(1, n), 1, 1:n);
                function rec(entries_sofar, i, remain)
                % Performs a recursion step
                %
                % Mutates the variables ``tableaux`` of the main function
                %
                % Args:
                %   entries_sofar (integer(1,\*)): Row vector of size ``n``, whose entries 1 to ``i`` have been populated
                %   i (integer): Current current entry to populate
                %   remain (integer(1,\*)): Numbers that can be filled in
                    if rowCount == self.dimension
                        return %stop when we know we've generated all tableaux
                    end
                    for j = remain % for loop over the numbers that we can put
                        if left_index(i) ==0 || entries_sofar(left_index(i)) < j
                            % either the i-th box doesn't have a left neighbor, or this neighbor index is less than the candidate j
                            if top_index(i)==0 || entries_sofar(top_index(i)) < j
                                % either the i-th box doesn't have a top neighbor, or this neighbor index is less than the candidate j
                                entries_sofar1 = entries_sofar;
                                entries_sofar1(i) = j;
                                if i == n
                                    rowCount = rowCount + 1;
                                    jPrimeFun(rowCount,entries_sofar1) = baseWords.word; 
                                    jRow(entries_sofar1) = baseWords.conjWord;
                                    self.basisOrder.(jRow) = rowCount;
                                    jFun(rowCount,:) = jRow;
                                else % i < n
                                    remain1 = setdiff(remain, j);
                                    rec(entries_sofar1, i + 1, remain1);
                                end
                            end
                        end
                    end
                end                
                function [above,left] = aboveLeft(part)
                % Calculates positional information of the
                % Young diagram
                %
                % Args:
                %  partition (integer(1,:)): Partition
                %
                % Returns:
                %  left (integer(1,:)): Index of the box immediately on the left in the Young diagram 0 if none
                %  above (integer(1,:)): Index of the box immediately to the top in the Young diagram, 0 if none
                    n = sum(part);
                    m = numel(part);
                    self.cSum = [0 cumsum(part(1:m-1))];
                    above = zeros(1,n);
                    left = 0:(n-1);
                    left(self.cSum+1)=0;
                    for j = 2:m
                        inds = 1:(part(j));
                        above(self.cSum(j)+inds) = self.cSum(j-1)+inds;
                    end
                end
            end
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
             rowFunEq = self.rowFunction(:,k) == self.rowFunction(:,k+1);  
             colFunEq = self.colFunction(:,k) == self.colFunction(:,k+1);  
             oneInds = self.rangeOfParts(rowFunEq);
             nOneInds = self.rangeOfParts(colFunEq);
             nInds = self.rangeOfParts(~rowFunEq&~colFunEq);
             len = numel(nInds);
             ratRows = zeros(1,2*len);
             ratCols = zeros(1,2*len);
             ratVals = zeros(1,2*len);
             count = 1;
             for i  = 1:len
                 if self.rowFunction(nInds(i),k)<self.rowFunction(nInds(i),k+1)
                    rPrime = self.rowFunction(nInds(i),transpotition(self,k));
                    rPrimeInd = self.basisOrder.(rPrime);
                    axDist = self.rowFunction(nInds(i),k+1)-self.rowFunction(nInds(i),k)+...
                    abs(self.colFunction(nInds(i),k+1)-self.colFunction(nInds(i),k));
                    ratRows(count:count+3) = [nInds(i) nInds(i) rPrimeInd rPrimeInd];
                    ratCols(count:count+3) = [nInds(i) rPrimeInd nInds(i) rPrimeInd];
                    ratVals(count:count+3) = [-1/axDist sqrt(1-1/axDist^2) sqrt(1-1/axDist^2) 1/axDist];
                    count = count + 4;
                 end
             end
            mat1 = sparse(ratRows,ratCols,ratVals,self.dimension,self.dimension);
            mat2 = sparse(oneInds,oneInds,1,self.dimension,self.dimension);
            mat3 = sparse(nOneInds,nOneInds,-1,self.dimension,self.dimension);
            im = mat1+mat2+mat3;
        end
        
        function rep = constructRep(self)
            n = self.group.domainSize;
            gens = cell(1,n-1);
            for i = 1:n-1
                gens{i} = self.transpotition(i);
            end
            self.group = self.group.subgroup(gens);
            images = cell(1,n-1);
            for i = 1:n-1
                images{i} = self.transImage(i);
            end
            rep = replab.RepByImages(self.group,self.field,self.dimension,images);
        end
        
        function t = transpotition(self,k)
            t = 1:self.group.domainSize;
            [t(k),t(k+1)] = deal(t(k+1),t(k));
        end
   end

end
