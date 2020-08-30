function [basis, basisMat] = findAllRepCoeffs2(rep,symb)
    n = rep.group.domainSize;
    [basisCell1,irreps] =  repBlockBasisEigAlg(rep,mult,symb);
    adjSwapIms = arrayfun(@(j) rep.image(replab.Permutation.transposition(n,j,j+1)),1:(n-1),'UniformOutput',false); 
    basis = arrayfun(@(i) extendBasis(basisCell1(i),irreps{i}.partition,adjSwapIms),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
end