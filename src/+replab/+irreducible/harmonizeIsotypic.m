function hi = harmonizeIsotypic(iso, context)
% Harmonizes an isotypic component
%
% As part of the operation, we identify the division algebra type (if the representations are over R),
% put those division algebras in their canonical basis, and make the irreducible
% representations not only equivalent but identical.
%
% Args:
%   iso (`+replab.Isotypic`): Isotypic component, not necessarily harmonized
%   context (`+replab.Context`): Sampling context
%
% Returns:
%   `+replab.HarmonizedIsotypic`: The harmonized isotypic component
    assert(isa(iso, 'replab.Isotypic'));
    assert(isa(context, 'replab.Context'));
    if isequal(iso.trivialDimension, iso.dimension)
        % trivial component, it's already harmonized
        hi = replab.HarmonizedIsotypic(iso.parent, iso.irreps, iso.E_internal);
        return
    end
    n = iso.nIrreps;
    if iso.overR
        irr1 = replab.irreducible.canonicalDivisionAlgebra(iso.irrep(1), context);
        W = irr1.A_internal;
        Winv = irr1.Ainv_internal;
    end
    irreps1 = cell(1, n);
    irreps1{1} = replab.rep.collapse(irr1);
    C = cell(1, n);
    C{1} = irr1.A_internal;
    for i = 2:n
        A = iso.changeOfBasis(1, i, context);
        Ainv = iso.changeOfBasis(i, 1, context);
        C{i} = W*A;
        irri = iso.irrep(i).similarRep(W*A, Ainv*Winv);
        irreps1{i} = replab.rep.collapse(irri);
    end
    E_internal = blkdiag(C{:}) * iso.E_internal;
    hi = replab.HarmonizedIsotypic(iso.parent, irreps1, E_internal);
end
