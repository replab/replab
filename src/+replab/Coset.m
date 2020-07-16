classdef Coset < replab.FiniteSet

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
    end

end
