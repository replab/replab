function sub = splitTrivialComponent(rep, context)
% Splits the trivial component out of a representation
    if isequal(rep.trivialDimension, 0) || isequal(rep.trivialDimension, rep.dimension)
        sub = replab.DispatchNext('Representation trivial component has already been split');
        return
    end
    T = rep.trivialSpace.sampleInContext(context, 1);
    % TODO: replace error estimate
    if ~replab.isNonZeroMatrix(T, replab.Parameters.doubleEigTol)
        s = replab.SubRep.fullSubRep(rep);
        rep.trivialDimension = 0;
        s.trivialDimension = 0;
        sub = {s};
        return
    end
    basis = orth(T);
    assert(size(basis, 1) == rep.dimension);
    dT = size(basis, 2);
    if rep.isUnitary
        [sub1 sub2] = rep.maschke(basis, basis');
    else
        [sub1 sub2] = rep.maschke(basis);
    end
    sub2.trivialDimension = 0;
    trivial = replab.Isotypic.fromTrivialSubRep(rep, sub1);
    sub = horzcat(trivial.irreps, {sub2});
end
