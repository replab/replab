function iso1 = refineIsotypic(iso, context)
% Implementation of `+replab.Isotypic.refine`
    if ~isequal(iso.isUnitary, true) || ~isequal(iso.parent.isUnitary, true)
        warning('Isotypic refinement only works for unitary representations');
        iso1 = iso;
        return
    end
    d = iso.irrepDimension;
    m = iso.multiplicity;
    P = full(iso.irrep(1).B_internal * iso.irrep(1).B_internal');
    for i = 2:m
        P = full(P + i*iso.irrep(i).B_internal * iso.irrep(i).B_internal');
    end
    P = iso.parent.commutant.project(P);
    [U D] = replab.numerical.sortedEig((P + P')/2, 'descend', false);
    start = 1;
    irreps = cell(1, m);
    for i = 1:m
        basis = U(:,start:start+d-1);
        irrep = iso.parent.subRep(basis, basis');
        irrep.isUnitary = true; % TODO: remove when dealing with nonunitary stuff
        irrep.isIrreducible = true;
        irrep.frobeniusSchurIndicator = iso.frobeniusSchurIndicator;
        irrep.isDivisionAlgebraCanonical = iso.isDivisionAlgebraCanonical;
        if ~isempty(iso.irrep(1).trivialDimension)
            irrep.trivialDimension = iso.irrep(1).trivialDimension;
        end
        irreps{i} = irrep;
        start = start + d;
    end
    A_internal = cell(1, m);
    Ainv_internal = cell(1, m);
    A_internal{1} = speye(d);
    Ainv_internal{1} = speye(d);
    C = iso.parent.commutant.sampleInContext(context, 1);
    for i = 2:m
        % A: space_i <- space_1
        A1 = irreps{i}.basis' * C * irreps{1}.basis;
        % We can use the other sample as well
        A2 = (irreps{1}.basis' * C * irreps{i}.basis)';
        A = A1 + A2';
        A = A * sqrt(d/real(trace(A*A')));
        A_internal{i} = A;
        Ainv_internal{i} = A';
    end
    iso1 = replab.Isotypic(iso.parent, irreps, A_internal, Ainv_internal);
end
