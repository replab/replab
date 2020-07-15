classdef CosetBase < replab.Obj

    properties (SetAccess = protected)
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
        subgroupChain % (`+replab.+bsgs.Chain`): Subgroup chain with base in lexicographic order
    end

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        subgroup % (`.FiniteGroup`): Subgroup of `.group`
    end

    methods

        function self = CosetBase(group, subgroup)
            assert(group.hasSameTypeAs(subgroup));
            self.group = group;
            self.subgroup = subgroup;
            self.isomorphism = group.niceMorphism;
            self.groupChain = self.isomorphism.imageGroup(group).lexChain;
            self.subgroupChain = self.isomorphism.imageGroup(subgroup).lexChain;
        end

    end

end
