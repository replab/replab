function irreps = identifyIrrepsInParent(sub, sample)
% Identifies the irreducible subrepresentations present in this subrepresentation
%
% If this subrepresentation should be first split further, it returns an empty array. Otherwise,
% it returns an array of one or two irreducible subrepresentations of `+replab.SubRep.parent`.
%
% Over the complex field (``sub.overC`` true), this merely checks that this subrepresentation is irreducible.
% If so, we set this subrepresentation ``sub.isIrreducible`` property to true, and we return it in
% a singleton cell array.
%
% Over the real field (``sub.overR`` true), three cases are possible:
%
% * This subrepresentation should be split further and is not irreducible in the sense below.
%   In that case, we return an empty cell array.
%
% * The subrepresentation does not encode a complex division algebra (``sub.divisionAlgebraName`` empty),
%   and we verify that this subrepresentation is irreducible. If so, we set this subrepresentation
%   ``sub.isIrreducible`` to true, and we return it in a singleton cell array.
%
% * The subrepresentation does encode a complex division algebra (``sub.divisionAlgebraName`` set to ``'C->R'``),
%   and moreover it splits into two complex conjugate subrepresentations. This method then identifies
%   the type of the irreducible subrepresentations (real-type, complex-type, quaternionic-type);
%   in the case of real-type subrepresentations, it splits them; in the case of quaternionic-type
%   subrepresentations, it identifies the canonical form of quaternion division algebra present.
%   It returns irreducible subrepresentations with ``sub.isIrreducible``, ``sub.divisionAlgebraName``
%   and ``sub.frobeniusSchurIndicator`` set to the relevant value.
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to identify
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%
% Returns:
%   cell(1,\*) of `+replab..SubRep`: Irreducible subrepresentations with their division algebra identified or ``[]``
    if sub.overR && strcmp(sub.divisionAlgebraName, 'C->R')
        if sub.isUnitary
            irreps = replab.irreducible.identifyIrrepsInParent_complexDivisionAlgebra_unitary(sub, sample);
        else
            irreps = replab.irreducible.identifyIrrepsInParent_complexDivisionAlgebra_nonunitary(sub, sample);
        end
    else
        irreps = replab.irreducible.identifyIrrepsInParent_trivialDivisionAlgebra(sub, sample);
    end
    if sub.inCache('trivialDimension') && sub.trivialDimension == 0
        for i = 1:length(irreps)
            irreps{i}.cache('trivialDimension', 0, '==');
        end
    end
end
