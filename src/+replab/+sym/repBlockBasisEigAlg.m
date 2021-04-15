function [bases,irreps] = repBlockBasisEigAlg(rep,symb)
    % Finds enough decomposition coefficients to block diagonalize a matrix in the commutant of a representation.
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
    %     basis (cell(1,*\) of double(*\,*\)): Bases of one multiplicity space,
    %         from each irrep
    %
    %     irreps: (cell(1,*\) of replab.sym.IntegerPartition):
    %         Partitions corresponding to the irreducible components of
    %         representation.

    n = rep.group.domainSize;
    parts = replab.sym.IntegerPartitions(n);
    mults = replab.ComplexCharacterTable.forPermutationGroup...
        (rep.group).multiplicities(rep);
    irrepInds = find(mults);
    irreps = parts.list(irrepInds);
    nIrreps = numel(irrepInds);
    csco = replab.sym.SymCSCO(n,0,rep.dimension,':');
    matList = csco.makeMatList(rep);
    eigenVals = csco.findSplitEigs(irreps);
    if ~symb
        if rep.isUnitary
            [V0,D0] = eig((matList{1}+matList{1}')/2,'vector');
        else
            [V0,D0] = eig(matList{1},'vector');
        end
        bases = arrayfun(@(index) getCoeffs(index),1:nIrreps,'UniformOutput',false);
    else
        stack = vertcat(matList{:});
        bases = arrayfun(@(index) getSymbCoeffs(index),1:nIrreps,'UniformOutput',false);
    end

    function cgs = getCoeffs(ind)
        splitEigs = eigenVals(ind,:);
        nEigs = numel(splitEigs);
        cgs =  V0(:,round(D0)==splitEigs(1));
        for i = 2:nEigs
            cgs = coeffFactors(i,splitEigs(i),cgs);
        end
    end

    function coeffs = getSymbCoeffs(ind)
        identityMatStack = sparse(1:csco.m*(rep.dimension),repmat(1:rep.dimension,1,csco.m),...
           repelem(eigenVals(ind,:),rep.dimension));
        coeffs = null(stack-identityMatStack,'r'); %replace with cycoltomic stuff here
    end

    function iso = coeffFactors(i,eigenvalue,oldIso)
        if rep.isUnitary
            mat = oldIso'*matList{i}*oldIso;
            [V,D] = eig((mat+mat')/2,'vector');
            reducedCG= V(:,round(D)==eigenvalue);
            iso = oldIso*reducedCG;
        else
            mat = oldIso\matList{i}*oldIso;
            [V,D] = eig(mat,'vector');
            reducedCG= V(:,round(D)==eigenvalue);
            iso = real(oldIso*reducedCG);
        end
    end
end