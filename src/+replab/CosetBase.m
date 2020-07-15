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
            if isa(group, 'replab.PermutationGroup')
                self.isomorphism = replab.FiniteIsomorphism.identity(group);
                self.groupChain = group.lexChain;
                self.subgroupChain = subgroup.lexChain;
            elseif isa(group, 'replab.NiceFiniteGroup')
                self.isomorphism = group.niceMorphism;
                self.groupChain = group.niceGroup.chain;
                self.subgroupChain = subgroup.niceGroup.chain;
            else
                error('unsupported');
            end
        end

    end

end
