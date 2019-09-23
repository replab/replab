function sub = decomposeExtractAllOnes(rep, E1, sampleE, sampleC)
% If a representation has the vector of all ones as the trivial component, we extract it
    if rep.isExtraFalse('hasTrivialSubspace')
        % trivial subspace has already been removed
        error('replab:dispatch:tryNext', 'try next');
    end 
    d = rep.dimension;
    field = rep.field;
    inSub = ones(d, 1);
    inParent = rep.U0' * rep.D0' * inSub;
    if ~E1.isEquivariant(inParent)
        % if the vector of all ones is not an invariant subspace, bail out
        error('replab:dispatch:tryNext', 'try next');
    end
    extra = struct;
    if rep.isExtraTrue('reducedBlocks')
        extra.reducedBlocks = true;
    end
    SB = replab.rep.standardBasis(d);
    UallOnes = rep.U0' * rep.D0' * SB(1,:) * rep.D0 * rep.U0;
    Urest = rep.U0' * rep.D0' * SB(2:end,:) * rep.D0 * rep.U0;
    extra1 = extra;
    extra1.hasTrivialSubspace = true;
    extra1.isIrreducible = true;
    if rep.overR
        extra1.divisionAlgebra = 'R';
    end
    allOnes = rep.subRepUnitary(UallOnes, extra1).collapseParent;
    rest = rep.subRepUnitary(Urest, extra).collapseParent;
    sub = {allOnes rest};
end
