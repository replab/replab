function sub = decomposeExtractTrivial(rep)
% If a representation has a trivial component, we extract it
    if rep.isExtraFalse('hasTrivialSubspace')
        % trivial subspace has already been removed
        error('replab:dispatch:tryNext', 'try next');
    end 
    d = rep.dimension;
    field = rep.field;
    extra = struct;
    if rep.isExtraTrue('reducedBlocks')
        extra.reducedBlocks = true;
    end
    U = replab.rep.standardBasis(d);
    % Extra for the trivial subspaces extracted
    extraT = extra;
    extraT.hasTrivialSubspace = true;
    extraT.isIrreducible = true;
    if rep.overR
        extraT.divisionAlgebra = 'R';
    end
    % Extra for the remaining subspaces
    extraR = extra;
    extraR.hasTrivialSubspace = false;
    % Use group averaging to find the trivial component
    Ed = rep.equivariant(rep.group.trivialRep(field, d));
    S = Ed.sample;
    if ~replab.isNonZeroMatrix(S, replab.Settings.doubleEigTol)
        sub = {rep.subRepUnitary(speye(d), extraR).collapseParent};
    else
        % Compute a basis of the trivial subspace
        Utrivial = orth(S)';
        % Compute the orthogonal subspace (nontrivial)
        Urest = null(Utrivial)';
        assert(size(Utrivial, 2) == d);
        assert(size(Urest, 2) == d);
        assert(size(Utrivial, 1) + size(Urest, 1) == d);
        % Construct the trivial irreducible subrepresentations
        dTrivial = size(Utrivial, 1);
        trivials = cell(1, dTrivial);
        for i = 1:dTrivial
            trivials{i} = rep.subRepUnitary(Utrivial(i,:), extraT).collapseParent;
        end
        % Construct the remaining representation
        rest = rep.subRepUnitary(Urest, extraR).collapseParent;
        sub = horzcat(trivials, {rest});
    end
end
