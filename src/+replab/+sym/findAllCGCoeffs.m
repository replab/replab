function [basis, basisMat, irreps] = findAllCGCoeffs(group,part1,part2,symb)
    % FInds all Clebsch Gordan Coefficients given of a tensor product
    % representation
    
    n = group.domainSize;
    rep = kron(group.irrep(part1,'orthogonal'),group.irrep(part2,'orthogonal'));
    [basisCell1,irreps] =  replab.sym.tensorBlockBasisEigAlg(rep,part1,part2,symb);
    adjSwapIms = arrayfun(@(j) sparse(rep.image(replab.Permutation.transposition(n,j,j+1))),1:(n-1),'UniformOutput',false); 
    basis = arrayfun(@(i) extendBasis(basisCell1(i),irreps{i}.partition,adjSwapIms),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
end