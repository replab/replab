function [basis, basisMat,irreps] = findAllRepCoeffs2(rep,isSymb)
    % Finds all decompositon coefficients given of a tensor product representation
    %
    % Args:
    %     group (replab.SymmetricGroup): Group being represented
    %     
    %     rep (replab.Rep): Representation to be decomposed. 
    %     
    %     symb (boolean): Do we want a symbolic result? Note that
    %         the result will be rational. Use seminormalToOrthogonal to
    %         help find the change of basis vectors to the orthogonal form.
    %     
    %     
    % Returns:
    %     basis (cell(1,*\) of cell(1,*\) of double(*\,*\)): Bases of all
    %         multiplicity spaces, organized by irrep, then by corresponding
    %         standard tableaux
    % 
    %     basisMat: (double(*\,*\)): change of basis matrix
    % 
    %     irreps: (cell(1,*\) of replab.sym.IntegerPartition):
    %         Partitions corresponding to the irreducible components of the
    %         representation
    
    n = rep.group.domainSize;
    [basisCell1,irreps] =  replab.sym.repBlockBasisEigAlg(rep,isSymb);
    numel(irreps)
    adjSwapIms = arrayfun(@(j) rep.image(replab.Permutation.transposition(n,j,j+1)),1:(n-1),'UniformOutput',false); 
    basis = arrayfun(@(i) replab.sym.extendBasis(basisCell1{i},irreps{i}.partition,adjSwapIms,isSymb),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
end