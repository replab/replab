function [bases,irreps] = tensorBlockBasis(rep,part1,part2)
    n = rep.group.domainSize;
    parts = replab.sym.findPartitions(n);
    mults = replab.CharacterTable.permutationCharTable...
        (rep.group).tensorProdMultiplicities(findPartIndices(part1,part2));
    irrepInds = find(flip(mults));
    irreps = parts.partCell(irrepInds);
    lat1 = replab.sym.YoungLattice(part1,n);
    s1 = lat1.symAndAntiSym;
    s2 = replab.sym.YoungLattice(part2,n).symAndAntiSym;
    subMatrixDim = s1(1)*s2(1)+s1(2)*s2(2);
    subMatInds = subMatrixIndices;
    nIrreps = numel(irrepInds);
    [xk,kGroups] = coeffs(n); %Coefficients for class operators C_2(k)
    m = numel(xk);
    stack = makeStack;
    bases = arrayfun(@(index) getCGCoeffs(index),1:nIrreps,'UniformOutput',false);
    
    function cgs = getCGCoeffs(ind)
        mult = mults(end-irrepInds(ind)+1);
        cgs = zeros(rep.dimension,mult);
        if irrepInds(ind) == 1
            cgs = oneDimCoeffs(false); %sign rep coefficents are known
            return
        elseif  irrepInds(ind) == parts.nParts
            cgs = oneDimCoeffs(true); %trival rep coefficents are also known
            return
        end
        identityMatStack = sparse(1:m*(subMatrixDim),repmat(1:subMatrixDim,1,m),...
           repelem(findSplitEigs(ind),subMatrixDim));
        subCG = nullForTest(stack-identityMatStack,mult);
        cgs(subMatInds,:) = subCG;
        function coeffs= oneDimCoeffs(isTrivial)
            if isTrivial
                dim = sum(s1);
                coeffs = zeros(dim^2,1);
                coeffs((dim*(0:(dim-1)))+(1:dim)) = lat1.generateTabFun;
            else
                [~,~,tabs] = lat1.generateTableaux;
                dim = sum(s1);
                signs = arrayfun(@(tableaux) replab.Permutation.sign(tabs(tableaux,:)),1:dim);
                coeffs = zeros(dim^2,1);
                coeffs((dim*(0:(dim-1)))+(dim:-1:1)) = signs;%orthToSemiBasis.*signs;
            end
        end
    end

    
    function Cstack = makeStack
        Cstack= zeros(m*(subMatrixDim),(subMatrixDim));
        for i = 1:m
            C0= zeros(subMatrixDim);
            for k = kGroups(i)+1:kGroups(i+1)
                C0 = C0 + xk{i}(k-kGroups(i))*D(2+k);
            end
            Cstack(((i-1)*(subMatrixDim)+1):(i*(subMatrixDim)),:)=C0;
        end
    end

    function eigenVals = findSplitEigs(ind)
         dEigs = dEigenvalues(n);
         eigenVals = zeros(1,m);
         for i = 1:m
             eigenVals(i)=xk{i}*dEigs(ind,kGroups(i)+1:kGroups(i+1))';
         end
    end

    function eigs = dEigenvalues(n)
        len =numel(irreps);
        eigs = zeros(len,n-2);
        for i = 1:len
            rec(parts.partCell{irrepInds(i)},n);
        end
        function rec(p,k)
            if sum(p) == 2
                return
            end
            row = p(end);
            col = numel(p);
            pSmall = p;
            pSmall(end) = pSmall(end)-1;
            if ~pSmall(end)
                pSmall = pSmall(1:(end-1));
            end
            eigs(i,k-2) = row-col;
            rec(pSmall,k-1);
        end
    end

   function [alpha,boundaries] = coeffs(n)
        %other choices will be given below each assignment for all n
        % which one to use depends on which is faster to have small
        % eigenvalues or smallar dimensions on the matrices
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
        len = numel(alpha);
        boundaries = zeros(1,len+1);
        boundaries(1)=0;
        for i = 1:len
            boundaries(i+1) = boundaries(i)+numel(alpha{i});
        end
    end

    function mat = D(k)
        mat = zeros(subMatrixDim);
        m2 = zeros(rep.dimension);
        for i = 1:(k-1)
            im = rep.image(swap(i,k));
            m2 = m2+ im;
            mat = mat + im(subMatInds,subMatInds);
        end
        function perm = swap(a,b)
            perm = 1:n;
            [perm(a),perm(b)] = deal(perm(b),perm(a));
        end
    end

    function inds = subMatrixIndices
        a = logical( [ ones(1,s1(1)) , zeros(1,s1(2)) ] );
        b = logical( [ ones(1,s2(1)) , zeros(1,s2(2)) ] );
        inds = kron(a,b)|kron(~a,~b);
    end

    

    function inds = findPartIndices(part1,part2)
        inds(1) = parts.nParts-parts.partitionHash.find(parts.toPowerForm(part1)')+1;
        inds(2) =  parts.nParts-parts.partitionHash.find(parts.toPowerForm(part2)')+1;
    end

    function kernel = nullForTest(mat,rank)
        [~,~,D] = svd(mat);
        kernel = D(:,(end-rank+1):end);
    end

end