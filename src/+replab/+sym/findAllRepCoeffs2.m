function [basis, basisMat] = findAllRepCoeffs2(rep,isRat)
    n = rep.group.domainSize;
    [basisCell1,irreps] =  replab.sym.repBlockBasisEigAlg(rep,isRat);
    adjSwapIms = arrayfun(@(j) rep.image(replab.Permutation.transposition(n,j,j+1)),1:(n-1),'UniformOutput',false); 
    basis = arrayfun(@(i) replab.sym.extendBasis(basisCell1(i),irreps{i}.partition,adjSwapIms),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
end