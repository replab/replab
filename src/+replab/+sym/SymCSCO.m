classdef SymCSCO
    % We borrowing the term CSCO from quantum mechanics.
    % Given a rep, this creates a commuting set of matrices with shared eigenvectors
    % Any 'basis vector' of the decomposition, (i.e. corresponding to a  
    % standard tableaux of an irreducible compoenent of the rep) is an
    % simultaneous eigenvector with a unique set of eigenvalues 
    %
    % This CSCO are the set of matrices A_k = sum_(i<k) image(swap(i,k))
    % or specially chosen linear combinations of these to reduce the number
    % of total commuting matrices
    %
    % See Section 4.6 of Chen's Group Representation Theory for Physicists
    properties
        forClebschGordan % bool: True if this is for clebsch-gordan decomposition
        domainSize %Integer: Domain size of symmetric group
        m % (integer): Number of matrices in the CSCO
    end
    
    properties(Access = protected)
        xk %cell of double(1,*\): How we linearly combine these matrices to 
        %make fewer total matrices in the CSCO
        % Each cell corresponds to a matrix in the CSCO and describes the
        % linear combination needed to get it.
        kGroups %double(1,*\): This is convinent for finding the A that each 
        % component in a cell corresponds to
        %Eg if xk was {[1 3] [-4 3 5] [4]}, kGroups would be [0 2 5 6]
        % So, we would know the first two A's are combined using the first cell
        % the 3rd, 4th and 5th are combined, as well, etc.
        %
        %For the next two, first know that, for CG coefficients, we only
        %need a sub matrix
        dim %(integer): Dimension of sub matrix, 
        % If not calculating CG,use  the dimension of the rep
        subMatIndices % (logical(*\)), Indices of sub matrix
        % If not calculating CG, use the character ':'
    end
    
    methods(Static)
        function eigs = aEigenvalues(n,integerPartitionList)
            len =numel(integerPartitionList);
            eigs = zeros(len,n-1);
            for i = 1:len
                rec(integerPartitionList{i}.partition,n);
            end

            function rec(p,k)
                if sum(p) == 1
                    return
                end
                row = p(end);
                col = numel(p);
                pSmall = p;
                pSmall(end) = pSmall(end)-1;
                if ~pSmall(end)
                    pSmall = pSmall(1:(end-1));
                end
                eigs(i,k-1) = row-col;
                rec(pSmall,k-1);
            end
        end
    end
    
    methods
        function self = SymCSCO(n,cg,dim,indices)
            self.domainSize = n;
            self.forClebschGordan = cg;
            [self.xk,self.kGroups] = self.coeffs;
            self.m = numel(self.xk);
            self.dim = dim;
            self.subMatIndices =indices; 
        end
        
        function eigenVals = findSplitEigs(self,partitionList)
             aEigs = self.aEigenvalues(self.domainSize,partitionList);
             eigenVals = zeros(numel(partitionList),self.m);
             for i = 1:self.m
                 eigenVals(:,i)=self.xk{i}*aEigs(:,self.forClebschGordan+(self.kGroups(i)+1:self.kGroups(i+1)))';
             end
        end

        function [alpha,boundaries] = coeffs(self)
            %other choices will be given below each assignment for all n
            % which one to use depends on which is faster to have small
            % eigenvalues or smallar dimensions on the matrices
            n = self.domainSize;
            if self.forClebschGordan
                switch n
                    case 3
                    alpha = {[1]};
                    case 4
                    alpha = {[0    -1]};
                    case 5
                    alpha = {[2    -3     1]};
                    case 6
                    alpha = {[-4    12     4     1]};
                    %{[0   1]  {3     1]};
                    %  {[2    -3     1], 1};
                    case 7
                    alpha ={[-16    -9    -5   -17    29]};
                    % {[2    -3     1], [-5    4]};
                    case 8
                    alpha = {[-7     2    -1], [36   -31    -4]};
                    %{[-4    12     4     1]  [7    -6]}
                    % {[-16    -9    -5   -17    29] 1};
                    case 9
                    alpha =  {[-16    -9    -5   -17    29], [8 -7]};
                        %gt{[-4    12     4     1]   [-60    54     5]}
                    % 
                    otherwise
                    % look at the best case for 8 and 9
                    % change the specifics depending on what the best choice is
                    if ~mod(n,2)
                        alpha = {[-4    12     4     1]  [7    -6]}; %best choice for 8 here
                        for i = 10:2:n
                            alpha{end+1} =[-i+1 i-2];
                        end
                    else
                        alpha ={[-16    -9    -5   -17    29], [8 -7]}; %best choice for 9 here
                        for i = 11:2:n
                            alpha{end+1} =[-i+1 i-2];
                        end
                    end   
                end
            else
                switch n
                    case 2
                        alpha = {1};
                    case 3
                    alpha = {[1    -1]};
                    case 4
                    alpha = {[4    0     -1]};
                    case 5
                    alpha = {[5    -5     2    -1]};
                    case 6
                    alpha = {[-15    16    -1     1     7]};
                    % {[4    0     -1] [-1 6]}
                    case 7
                    alpha ={ [-15    16    -1     1     7], 1};
                    % {[10    49     9   -12     3   -86]}
                    % {[5    -5     2    -1], [7 -6]}
                    %
                    case 8
                    alpha =  {[-15    16    -1     1     7],  [-8     7]};
                    %{[-15    16    -1     1     7], 1 1}
                    case 9
                    alpha = {[-15    16    -1     1     7], [50   -44    -5]};
                    % {[10    49     9   -12     3   -86] [-9     8]}
                    % {[5    -5     2    -1], [7 -6], [9 -8]}
                    otherwise
                    % look at the best case for 8 and 9
                    % change the specifics depending on what the best choice is
                    if mod(n,2)
                        alpha = {[-15    16    -1     1     7],  [-8     7]}; %best choice for 8 here
                        for i = 10:2:n
                            alpha{end+1} =[-i i-1];
                        end
                    else
                        alpha = {[-15    16    -1     1     7], [50   -44    -5]}; %best choice for 9 here
                        for i = 11:2:n
                            alpha{end+1} =[-i i-1];
                        end
                    end   
                end
            end
            len = numel(alpha);
            boundaries = zeros(1,len+1);
            boundaries(1)=0;
            for i = 1:len
                boundaries(i+1) = boundaries(i)+numel(alpha{i});
            end
        end

        function list = makeMatList(self,rep)
            list=cell(1,self.m);
            zerothK = 1 +self.forClebschGordan;
            for i = 1:self.m
                C0= zeros(self.dim);
                for k = self.kGroups(i)+1:self.kGroups(i+1)
                    C0 = C0 + self.xk{i}(k-self.kGroups(i))*self.A(zerothK+k,rep);
                end
                list{i}=C0;
            end
        end
        
         function mat = A(self,k,rep)
            mat = zeros(self.dim);
            for i = 1:(k-1)
                im = rep.image(replab.Permutation.transposition(self.domainSize,i,k));
                mat = mat +im(self.subMatIndices,self.subMatIndices);
            end
        end
    end
end