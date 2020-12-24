function hi = Isotypic_harmonize(iso, context)
% Harmonizes an isotypic component
%
% As part of the operation, we identify the division algebra type (if the representations are over R),
% put those division algebras in their canonical basis, and make the irreducible
% representations not only equivalent but identical.
%
% Args:
%   iso (`+replab.Isotypic`): Isotypic component with ``isHarmonized`` false and ``trivialDimension == 0``
%   context (`+replab.Context`): Sampling context
%
% Returns:
%   `+replab.Isotypic`: The harmonized isotypic component
    assert(isa(iso, 'replab.Isotypic'));
    assert(~iso.isHarmonized && iso.trivialDimension == 0);
    n = iso.nIrreps;
    irreps = cell(1, n);
    P = iso.projection_internal;
    if iso.overR
        irr1 = replab.irreducible.canonicalDivisionAlgebra(iso.irrep(1), context);
        W = irr1.A_internal;
        Winv = irr1.Ainv_internal;
        irreps{1} = irr1.collapse;
        range = iso.irrepRange(1);
        P(range, :) = W * P(range, :);
    else
        W = speye(iso.irrepDimension);
        Winv = W;
        irreps{1} = iso.irrep(1);
    end
    for i = 2:n
        [A Ainv] = replab.irreducible.Isotypic_changeOfBasis(iso, 1, i, context);
        irri = iso.irrep(i).similarRep(W * A, 'inverse', Ainv * Winv).collapse;
        range = iso.irrepRange(i);
        P(range, :) = W * A * P(range, :);
        if irreps{1}.inCache('isUnitary')
            irri.cache('isUnitary', irreps{1}.isUnitary, '==');
        end
        irri.cache('frobeniusSchurIndicator', irreps{1}.frobeniusSchurIndicator, '==');
        irri.cache('isDivisionAlgebraCanonical', irreps{1}.isDivisionAlgebraCanonical, '==');
        irreps{i} = irri;
    end
    isHarmonized = true;
    hi = replab.Isotypic(iso.parent, irreps, P, iso.irrepDimension, isHarmonized);
end
