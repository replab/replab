function irreps = decomposeExtractAllOnes(rep)
% If a representation has the vector of all ones as the trivial component, we extract it
    if rep.isExtraFalse('hasTrivialSubspace')
        % trivial subspace has already been removed
        error('replab:dispatch:tryNext', 'try next');
    end 
    d = rep.dimension;
    field = rep.field;
    E1 = rep.equivariant(rep.group.trivialRep(field, 1));
    if ~E1.isEquivariant(ones(d, 1))
        % if the vector of all ones is not an invariant subspace, bail out
        error('replab:dispatch:tryNext', 'try next');
    end
    extra = struct;
    if rep.isExtraTrue('reducedBlocks')
        extra.reducedBlocks = true;
    end
    U = replab.rep.standardBasis(d);
    extra1 = extra;
    extra1.hasTrivialSubspace = true;
    extra1.isIrreducible = true;
    if rep.overR
        extra1.divisionAlgebra = 'R';
    end
    irrep1 = rep.subRepUnitary(U(1,:), extra1).collapseParent;
    rest = rep.subRepUnitary(U(2:end,:), extra).collapseParent;
    irreps = horzcat({irrep1}, replab.rep.decompose(rest));
end
