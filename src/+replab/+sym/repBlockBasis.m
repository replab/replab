function [bases,irreps] = repBlockBasis(rep)
    n = rep.group.domainSize;
    parts = replab.sym.findPartitions(n);
    mults = replab.CharacterTable.permutationCharTable...
        (rep.group).multiplicities(rep);
    irrepInds = find(flip(mults));
    irreps = parts.partCell(irrepInds);
    nIrreps = numel(irrepInds);
    [xk,kGroups] = coeffs(n); %Coefficients for class operators C_2(k)
    m = numel(xk);
    stack = makeStack;
    bases = arrayfun(@(index) getCoeffs(index),1:nIrreps,'UniformOutput',false);
    
    function coeffs = getCoeffs(ind)
        mult = mults(end-irrepInds(ind)+1);
        identityMatStack = sparse(1:m*(rep.dimension),repmat(1:rep.dimension,1,m),...
           repelem(findSplitEigs(ind),rep.dimension));
        coeffs = nullForTest(stack-identityMatStack,mult);
    end

    function Cstack = makeStack
        Cstack= zeros(m*(rep.dimension),rep.dimension);
        for i = 1:m
            C0= zeros(rep.dimension);
            for k = kGroups(i)+1:kGroups(i+1)
                C0 = C0 + xk{i}(k-kGroups(i))*D(1+k);
            end
            Cstack(((i-1)*rep.dimension+1):(i*rep.dimension),:)=C0;
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
        eigs = zeros(len,n-1);
        for i = 1:len
            rec(parts.partCell{irrepInds(i)},n);
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

    function [alpha,boundaries] = coeffs(n)
        %other choices will be given below each assignment for all n
        % which one to use depends on which is faster to have small
        % eigenvalues or smallar dimensions on the matrices
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
            alpha = {[-15    16    -1     1     7], 1 1};
            % {[-15    16    -1     1     7],  [-8     7]}
            case 9
            alpha = {[10    49     9   -12     3   -86] [-9     8]};
            % {[-15    16    -1     1     7], [50   -44    -5]}
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
        len = numel(alpha);
        boundaries = zeros(1,len+1);
        boundaries(1)=0;
        for i = 1:len
            boundaries(i+1) = boundaries(i)+numel(alpha{i});
        end
    end

    function mat = D(k)
        mat = zeros(rep.dimension);
        for i = 1:(k-1)
            mat = mat + rep.image(swap(i,k));
        end
        function perm = swap(a,b)
            perm = 1:n;
            perm([a b]) = perm([b a]);
        end
    end

    function kernel = nullForTest(mat,rank)
        [~,~,D] = svd(mat);
        kernel = D(:,(end-rank+1):end);
    end

end