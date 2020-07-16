classdef Coset < replab.FiniteSet

    properties (SetAccess = protected)
        parent % (`.FiniteGroup`): Group containing this coset
        group % (`.FiniteGroup`): Group
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
    end

end
