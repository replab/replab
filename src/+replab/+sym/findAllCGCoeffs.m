function [basis, basisMat, irreps] = findAllCGCoeffs(group,part1,part2,isSymb)
    % FInds all Clebsch Gordan coefficients of a tensor product
    % representation
    % Args:
        % group (replab.SymmetricGroup): Group being represented
        %
        % part1 (integer(*\)): Partition of first irrep in tensor product
        %
        % part2 (integer(*\)): Partition of second irrep in tensor product
        %
        % isSymb (boolean): Do we want a symbolic result? Note that
        % the result will be rational. Use seminormalToOrthogonal to
        % help find the change of basis vectors to the orthogonal form.
        %
        %
    % Returns:
        % basis (cell(1,*\) of cell(1,*\) of double(*\,*\)): Bases of all
        % multiplicity spaces, organized by irrep, then by corresponding standard
        % tableaux
        %
        % basisMat: (double(*\,*\)): change of basis matrix
        %
        % irreps: (cell(1,*\) of replab.sym.IntegerPartition):
        % Partitions corresponding to the irreducible components of the
        % tensor product
    n = group.domainSize;
    if ~isSymb
        rep = kron(group.irrep(part1,'orthogonal'),group.irrep(part2,'orthogonal'));
    else
        rep = kron(group.irrep(part1,'seminormal'),group.irrep(part2,'seminormal'));
    end
    [basisCell1,irreps] =  replab.sym.tensorBlockBasisEigAlg(rep,part1,part2,isSymb);
    adjSwapIms = arrayfun(@(j) sparse(rep.image(replab.Permutation.transposition(n,j,j+1))),1:(n-1),'UniformOutput',false); 
    %We do this because we know the ajacent transpotions are sparse
    basis = arrayfun(@(i) replab.sym.extendBasis(basisCell1{i},irreps{i}.partition,adjSwapIms,...
        isSymb),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
end