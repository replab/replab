classdef Cosets < replab.Obj

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        subgroup % (`.FiniteGroup`): Subgroup of `.group`
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
        subgroupChain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
    end

    methods

        function self = Cosets(group, subgroup)
            assert(group.hasSameTypeAs(subgroup));
            self.group = group;
            self.subgroup = subgroup;
            self.isomorphism = group.niceMorphism;
            self.groupChain = group.niceGroup.lexChain;
            self.subgroupChain = subgroup.niceGroup.lexChain;
        end

    end

end
