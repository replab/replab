classdef SymmetricSpechtIrrep < replab.Rep
% An irreducible representation of a symmetric group
%
% Each irrep corresponds to an unordered partition of an integer ``n``.
%
% Ex: ``4 = 2+2`` and ``4 = 3+1`` both generate an irrep of $S_4$.
%
% The entries of images in these represntations consist of only ``0``, ``1``, and ``-1``.
%
% This implementation is based on
% Wiltshire-Gordon, John D.; Woo, Alexander; Zajaczkowska, Magdalena (2017).
% "Specht Polytopes and Specht Matroids", `<https://arxiv.org/abs/1701.05277>`_

    properties (SetAccess = protected)
       partition % (integer(1,\*)): The generating partition of ``n``
       conjugatePartition % (integer(1,\*)): The conjugate partition is the partition obtained by transposing the Young diagram of a partition
       indepRowWords % (integer(\*,\*)): The row words; corresponding to all standard tableaux
       indepColWords % (integer(\*,\*)): The column words; corresponding to all standard tableaux
       basis % (integer(\*,\*)): The submatrix of the Specht matrix whose rows and columns correspond to the row and column words of standard tableaux
    end

    properties (Access=protected)
        cSum
        underlyingRep % This is an underlying RepByImages object used to quickly find the image
    end

   methods
       function self = SymmetricSpechtIrrep(group, part)
        % Constructs an irreducible representation of S_n
        %
        % Args:
        %   partition (integer(1,\*)): Partition
        %   group (`+replab.Group`): Symmetric group being represented
        %
        % Returns:
       %    `+replab.Rep`: Specht representation of the group
            if sum(part) ~= group.domainSize
                error('This is not a valid Young Diagram for this domain size.')
            end
            self.field = 'R';
            self.group = group;
            self.partition = part;
            self.dimension = replab.sym.partition.dimension(part);
            self.conjugatePartition = replab.sym.partition.conjugatePart(part);
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

   methods (Access = protected)

       function symIrrepHelper(self)
       % Helper function for constructor
           len = numel(self.partition);
           self.cSum = [0 cumsum(self.partition(1:len-1))];
           baseWords = replab.sym.words(self.partition,self.conjugatePartition);
           %Generate our choice of base words corresponding to the partition and
           %conjuagate partition.
           [self.indepRowWords,self.indepColWords] = tabs(self.partition,...
                                                          self.dimension,baseWords.word,baseWords.conjWord);
           %Generate words corresponding to linearly independent
           %columns and rows
           self.basis = self.subSpecht(self.indepRowWords);
           function [rowWords,colWords] = tabs(partition,dim,word,cWord)
           % Enumerate standard Young tableaux for a given partition
           %
           % Indices in the Young diagram increase left to right then top to bottom
           %
           % These tableaux are in one-to-one correspondence with rearrangements of words
           %
           % Returns:
           %  rowWords (integer(:,:): A matrix whose rows enumerate all row words corresponding to standard Young tableaux
           %  colWords (integer(:,:): A matrix whose rows enumerate all column words corresponding to standard Young tableaux
               n = sum(partition);
               rowWords = zeros(dim,n);
               colWords = zeros(dim,n);
               rowCount = 0;
               [left_index, top_index] = aboveLeft(partition);
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
                   for j = remain % for loop over the numbers that we can put
                       if left_index(i) ==0 || entries_sofar(left_index(i)) < j
                           % either the i-th box doesn't have a left neighbor, or this neighbor index is less than the candidate j
                           if top_index(i)==0 || entries_sofar(top_index(i)) < j
                               % either the i-th box doesn't have a top neighbor, or this neighbor index is less than the candidate j
                               entries_sofar1 = entries_sofar;
                               entries_sofar1(i) = j;
                               if i == n
                                   rowCount = rowCount + 1;
                                   colWords(rowCount,entries_sofar1) = word;
                                   rowWords(rowCount,entries_sofar1) = cWord;
                               else % i < n
                                   remain1 = setdiff(remain, j);
                                   rec(entries_sofar1, i + 1, remain1);
                               end
                           end
                       end
                   end
               end

               function [above,left] = aboveLeft(part)
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
               perm = self.cSum(col) + row;
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