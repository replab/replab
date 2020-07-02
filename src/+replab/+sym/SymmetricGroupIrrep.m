classdef SymmetricGroupIrrep < replab.Rep
    
    properties (SetAccess = protected)
        partition
        conjugatePartition
        isEqualToConjugateRep
    end


    properties (GetAccess = protected, SetAccess = protected)        
        nDiagRows
        nDiagCols
        word
        cWord
        dimWord
        rowHashTable
        colHashTable
        specht
        rowWordList
        underlyingRepObject
     end
    
     methods(Static)
        function conj = conjugatePart(part)
            m = max(part);
            conj = zeros(1,m);
            for j = 1:m
                conj(j) = sum(j-1<part);
            end
        end
        
        function dim = dimension(part)
            n = sum(part);
            [~, cWord,dimWord] = replab.sym.SymmetricGroupIrrep.words(part,replab.sym.SymmetricGroupIrrep.conjugatePart(part));
            columns = zeros(1,n);
            for k = 1:n
                columns(k) = sum(cWord(k+1:n)==cWord(k));
            end    
            dim = round(factorial(n)/prod(dimWord+columns));
        end     
     end
     
    methods(Access = protected, Static)
        function [word, cWord,dimWord] = words(part,cPart)
            N = length(part);
            l = sum(part);
            word = repelem(1:N,part);
            cWord = zeros(1,l);
            dimWord = zeros(1,l);
            index = 1;
            for m = 1:max(cPart)
                p = sum(cPart > m - 1);
                cWord(index:index+p-1) = 1:p;
                dimWord(index:index+p-1) = p:-1:1;
                index = index + p;
            end
        end
        
        function h = hashFun(word,m,n)
            h = (word-1)*(m.^(0:(n-1)))'+1;
        end
     end

     methods
         function self = SymmetricGroupIrrep(group, part)
             if group.domainSize > 12
                 error('The domain size must be less than or equal to 12.')
             end
             if group.domainSize > 8
                 warning('If the domain size exceeds 8, construction may be slow. For n > 9, some partitions may not work.')
             end
             if sum(part) ~= group.domainSize
                 error('This is not a valid Young Diagram for this domain size.')
             end
             self.group = group;
             self.field = 'R';
             self.isIrreducible = true;
             self.nDiagCols = max(part);
             self.nDiagRows = length(part);
             self.partition = part;
             self.conjugatePartition = self.conjugatePart(part);
             self.isEqualToConjugateRep = isequal(self.partition,self.conjugatePartition);
             [self.word, self.cWord, self.dimWord] = self.words(part,self.conjugatePartition);
             %%Dimension Calculation
             n = group.domainSize;
             columns = zeros(1,n);
             for k = 1:n
                columns(k) = sum(self.cWord(k+1:n)==self.cWord(k));
             end
             self.dimension = factorial(n)/prod(self.dimWord+columns);
             if self.nDiagRows == 1
                 self.underlyingRepObject = group.trivialRep('R',1);
                 self.isUnitary = true;
                 self.trivialDimension = 1;
             elseif self.nDiagCols == 1
                 self.underlyingRepObject = group.signRep();
                 self.isUnitary = true;
                 self.trivialDimension = 0;
             elseif self.nDiagRows == 2 && self.nDiagCols == n-1
                 self.underlyingRepObject = group.standardRep();
                 self.trivialDimension = 0;
             elseif self.nDiagCols == 2 && self.nDiagRows == n-1
                 self.underlyingRepObject = group.tensorRep(self.field,{group.standardRep(),group.signRep()});
                 self.trivialDimension = 0;
             else   
                 nRows = factorial(n)/prod(factorial(self.conjugatePartition));
                %nCols = factorial(n)/prod(factorial(part));
                try
                    self.rowWordList = zeros(nRows,n);
                    self.rowHashTable = zeros(1,self.nDiagCols^n);
                    self.colHashTable = zeros(1,self.nDiagRows^n);
                catch
                    error('The partition does not work in this implementation.')
                end
                
                genPerms;
                self.underlyingRepObject = constructRBIObject;
                self.trivialDimension = 0;
             end

        
            function genPerms
            spechtEntries = zeros(1,factorial(n));
            spechtRows = zeros(1,factorial(n));
            spechtCols = zeros(1,factorial(n));
            arr = 1:n; 
            sign = 1;
            permCount = 1;
            newSpechtCol = 1;
            newSpechtRow = 1; 
            rec(n);
            function rec(k)
            parity = mod(k,2);
            if k == 1
                updateSpecht();
                permCount = permCount + 1;
            else
                rec(k-1);
                for i = 1:(k-1)
                    if parity == 0
                        arr([i k]) = arr([k i]);
                        sign = sign *-1;
                    else
                        arr([1 k]) = arr([k 1]);
                        sign = sign *-1;
                    end
                rec(k-1);
                end
            end
            function updateSpecht()
                cWord1 = self.cWord(arr);
                word1 = self.word(arr);
                hash = self.hashFun(word1,self.nDiagRows,n);
                cHash = self.hashFun(cWord1,self.nDiagCols,n);
                if self.rowHashTable(cHash)==0 && self.colHashTable(hash)==0
                    spechtRows(permCount)  = newSpechtRow;
                    spechtCols(permCount)  = newSpechtCol;
                    spechtEntries(permCount)  = sign;
                    self.specht(newSpechtRow,newSpechtCol) = sign;
                    self.rowHashTable(cHash) = newSpechtRow;
                    self.rowWordList(newSpechtRow,:) = cWord1;
                    self.colHashTable(hash) = newSpechtCol;
                    newSpechtRow = newSpechtRow + 1;
                    newSpechtCol = newSpechtCol +1;
                elseif self.rowHashTable(cHash)==0
                    spechtRows(permCount)  = newSpechtRow;
                    spechtCols(permCount)  = self.colHashTable(hash);
                    spechtEntries(permCount)  = sign;
                    self.rowHashTable(cHash) = newSpechtRow;
                    self.rowWordList(newSpechtRow,:) = cWord1;
                    newSpechtRow = newSpechtRow + 1;
                elseif self.colHashTable(hash)==0
                    spechtRows(permCount)  = self.rowHashTable(cHash);
                    spechtCols(permCount)  = newSpechtCol;
                    spechtEntries(permCount)  = sign;
                    self.colHashTable(hash) = newSpechtCol;
                    newSpechtCol = newSpechtCol + 1;
                else
                    spechtRows(permCount)  = self.rowHashTable(cHash);
                    spechtCols(permCount)  = self.colHashTable(hash);
                    spechtEntries(permCount)  = sign;
                end
             end
        
            end
            self.specht = sparse(spechtRows,spechtCols,spechtEntries);
            end
             
            function irrep = constructRBIObject
                    %This speedup to finding pivot columns was found on
                    %www.mathworks.com/matlabcentral/fileexchange/21583-fast-reduced-row-echelon-form.
                    R = qr(self.specht);
                    [indep_rows, colPiv] = max(R ~= 0, [], 2); 
                    indep_rows = full(indep_rows); 
                    colPiv = colPiv(indep_rows);
                    basis = self.specht(:,colPiv);
                    irrep = replab.RepByImages.fromImageFunction(self.group, 'R', self.dimension, @(g) naiveImage(g));
                function imag = naiveImage(perm)
                     cPerm = self.rowHashTable(self.hashFun(self.rowWordList(:,perm),self.nDiagCols,n));
                     action = self.specht(cPerm,colPiv);
                     imag = round(basis\action);
                end
             end

        end
         
        function rho = image_internal(self, g)
            rho = self.underlyingRepObject.image(g);
        end
        
     end
end