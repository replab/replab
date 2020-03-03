function hi = harmonizeIsotypic(iso, context)
% Harmonizes an isotypic component
    assert(isa(iso, 'replab.Isotypic'));
    assert(isa(context, 'replab.Context'));
    if isequal(iso.trivialDimension, iso.dimension)
        % trivial component, it's already harmonized
        hi = replab.HarmonizedIsotypic(iso.parent, iso.irreps);
        return
    end
    n = iso.nIrreps;
    C = iso.parent.commutant.sampleInContext(context, 1);
    if iso.overR
        irr1 = replab.irreducible.canonicalDivisionAlgebra(iso.irrep(1), context);
        W = irr1.A_internal;
        Winv = irr1.Ainv_internal;
    end
    irreps1 = cell(1, n);
    irreps1{1} = replab.rep.collapse(irr1);
    for i = 2:n
        A = iso.changeOfBasis(1, i, context);
        Ainv = iso.changeOfBasis(i, 1, context);
        irri = iso.irrep(i).similarRep(W*A, Ainv*Winv);
        irreps1{i} = replab.rep.collapse(irri);
    end
    hi = replab.HarmonizedIsotypic(iso.parent, irreps1);
end
