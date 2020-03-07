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
        assert(isequal(irr1.isUnitary, true));
        iso1 = iso.changeIrrepBasis(1, W, Winv);
    else
        iso1 = iso;
    end
    A_list = cell(1, n);
    Ainv_list = cell(1, n);
    for i = 1:n
        [A Ainv] = iso1.changeOfBasis(1, i, context);
        Ainv = inv(A);
        A_list{i} = A;
        Ainv_list{i} = Ainv;
    end
    iso2 = iso1.changeEachIrrepBasis(A_list, Ainv_list);
    for i = 2:n
        irrepsi = iso2.irreps{i};
        irrepsi.isUnitary = iso2.irreps{1}.isUnitary;
        irrepsi.frobeniusSchurIndicator = iso2.irreps{1}.frobeniusSchurIndicator;
        irrepsi.isDivisionAlgebraCanonical = iso2.irreps{1}.isDivisionAlgebraCanonical;
    end
    hi = replab.HarmonizedIsotypic(iso2.parent, iso2.irreps, iso2.E_internal);
end
