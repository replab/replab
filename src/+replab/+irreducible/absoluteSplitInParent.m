function subs = absoluteSplitInParent(sub, sample, forceNonUnitaryAlgorithms)
% Decomposes this subrepresentation into subrepresentations of its parent
%
% It tries to identify irreducible subrepresentations but does not certify their irreducibility.
%
% Returned subrepresentations are subrepresentations of ``sub.parent``, so the process can be iterated if necessary.
%
% In the complex case (``sub.overC`` true), it returns subrepresentations that should be irreducible. The irreducibility
% of those subrepresentations should then be certified using a relevant algorithm.
%
% In the real case (``sub.overR`` true), it returns subrepresentations that fall in one of these two cases:
%
% * The subrepresentation encodes a pair of supposedly irreducible complex subrepresentations. In that case,
%   ``subs{i}.divisionAlgebraName`` is set to ``'complex'``, and the subrepresentation should be processed further
%   to identify the type of the real irreducible representation(s) contained.
%
% * The subrepresentation should be irreducible and of real-type. In that case, ``subs{i}.divisionAlgebraName`` is
%   empty. The subrepresentation should directly be certified irreducible using a relevant algorithm.
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to split
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%   forceNonUnitaryAlgorithms (logical): Whether to force the use of algorithms for not necessarily unitary representations
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of this subrepresentation
    if sub.dimension <= 1
        subs = {sub};
        return
    end
    if sub.isUnitary && ~forceNonUnitaryAlgorithms
        if sub.overC
            subs = replab.irreducible.absoluteSplitInParent_complex_unitary(sub, sample);
        else
            subs = replab.irreducible.absoluteSplitInParent_real_unitary(sub, sample);
        end
    else
        if sub.overC
            subs = replab.irreducible.absoluteSplitInParent_complex_nonunitary(sub, sample);
        else
            subs = replab.irreducible.absoluteSplitInParent_real_nonunitary(sub, sample);
        end
    end
    if sub.inCache('trivialDimension') && sub.trivialDimension == 0
        for i = 1:length(subs)
            subs{i}.cache('trivialDimension', 0, '==');
        end
    end
end
