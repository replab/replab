function [A Ainv] = Isotypic_changeOfBasis(iso, i, j, context)
% Returns change of basis matrices that relate two irreducible representations
%
% ``A`` such that ``A * iso.irrep(j).image(g) * Ainv = iso.irrep(i).image(g)``
%
% Args:
%   i (integer): Index of an irreducible representation
%   j (integer): Index of an irreducible representation
%   context (`+replab.Context`): Sampling context
%
% Returns
% -------
%   A: double(\*,\*), may be sparse
%     Change of basis matrix
%   Ainv: double(\*,\*), may be sparse
%     Inverse of change of basis matrix
    if i == j
        A = speye(iso.irrepDimension);
        Ainv = A;
        return
    end
    C = iso.parent.commutant.sampleInContext(context, 1);
    A = iso.irrep(i).projection_internal * C * iso.irrep(j).injection_internal;
    A = A * sqrt(iso.irrepDimension/real(trace(A*A'))) * sign(A(1,1));
    if iso.irrep(i).isUnitary && iso.irrep(j).isUnitary
        Ainv = A';
    elseif iso.overC || isequal(iso.irrep(1).frobeniusSchurIndicator, 1)
        Ainv = iso.irrep(j).projection_internal * C * iso.irrep(i).injection_internal;
        Ainv = Ainv/(trace(A*Ainv)/iso.irrepDimension);
    else
        Ainv = inv(A);
    end
end
