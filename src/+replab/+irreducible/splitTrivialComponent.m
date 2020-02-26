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
    sub = cell(1, dT+1);
    B = sub1.B_internal;
    E = sub1.E_internal;
    for i = 1:dT
        sub{i} = rep.subRep(B(:,i), E(i,:));
        sub{i}.isIrreducible = true;
        sub{i}.trivialDimension = 1;
        sub{i}.frobeniusSchurIndicator = 1;
    end
    sub{dT+1} = sub2;
end
