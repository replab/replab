function [bases,irreps] = tensorBlockBasisEigAlg(rep,part1,part2,symb)
    n = rep.group.domainSize;
    parts = replab.sym.IntegerPartitions(n);
    mults =  replab.CharacterTable.permutationCharTable(rep.group).tensorProdMultiplicities...
        ([parts.index(part1) parts.index(part2)]);
    irrepInds = find(mults);
    irreps = parts.list(irrepInds);
    lat1 = replab.sym.YoungLattice(part1,n);
    s1 = lat1.symAndAntiSym;
    s2 = replab.sym.YoungLattice(part2,n).symAndAntiSym;
    subMatrixDim = s1(1)*s2(1)+s1(2)*s2(2);
    subMatInds = subMatrixIndices;
    nIrreps = numel(irrepInds);
    csco = symCSCO(n,1,subMatrixDim,subMatInds);
    matList = csco.makeMatList;
    eigenVals = csco.findSplitEigs(irreps);
    if ~symb
        if rep.isUnitary
            [V0,D0] = eig((matList{1}+matList{1}')/2,'vector');
        else
            [V0,D0] = eig(matList{1},'vector');
        end
        bases = arrayfun(@(index) getCGCoeffs(index),1:nIrreps,'UniformOutput',false);
    else
        stack = vertcat(matList{:});
        bases = arrayfun(@(index) getSymbCGCoeffs(index),1:nIrreps,'UniformOutput',false);
    end
    
    
    function cgs = getCGCoeffs(ind)
        mult = mults(irrepInds(ind));
        cgs = zeros(rep.dimension,mult);
        if irrepInds(ind) == parts.nParts
            cgs = oneDimCoeffs(false); %sign rep coefficents are known
            return
        elseif  irrepInds(ind) == 1
            cgs = oneDimCoeffs(true); %trival rep coefficents are also known
            return
        end
        splitEigs = eigenVals(ind,:);
        nEigs = numel(splitEigs);
        subCG =  V0(:,round(D0)==splitEigs(end));
        for i = 1:nEigs-1
            subCG = reducedClebschGordan(i,splitEigs(i),subCG);
        end
        cgs(subMatInds,:) = subCG;
        function coeffs= oneDimCoeffs(isTrivial)
            if isTrivial
                dim = sum(s1);
                coeffs = zeros(dim^2,1);
                coeffs((dim*(0:(dim-1)))+(1:dim)) = 1/sqrt(dim);
            else
                [~,~,tabs] = lat1.generateTableaux;
                dim = sum(s1);
                signs = arrayfun(@(tableaux) replab.Permutation.sign(tabs(tableaux,:)),1:dim);
                coeffs = zeros(dim^2,1);
                coeffs((dim*(0:(dim-1)))+(dim:-1:1)) = signs/sqrt(dim);%orthToSemiBasis.*signs;
            end
        end
    end

    function cgs = getSymbCGCoeffs(ind)
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
            subCG = nullForTest(stack-identityMatStack,mult); %You would do cyclotomic stuff here
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

    function iso = reducedClebschGordan(i,eigenvalue,oldIso)
        mat = oldIso'*matrixList{i}*oldIso;
        [V,D] = eig((mat+mat')/2,'vector');
        reducedCG= V(:,round(D)==eigenvalue);
        iso = oldIso*reducedCG;
    end
    
    function inds = subMatrixIndices
        a = logical( [ ones(1,s1(1)) , zeros(1,s1(2)) ] );
        b = logical( [ ones(1,s2(1)) , zeros(1,s2(2)) ] );
        inds = kron(a,b)|kron(~a,~b);
    end

end