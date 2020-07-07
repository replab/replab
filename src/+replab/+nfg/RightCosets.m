classdef RightCosets < replab.Str

    properties (SetAccess = protected)
        group % (`.NiceFiniteGroup`): Group
        subgroup % (`.NiceFiniteGroup`): Subgroup of `.group`
    end

    properties
        nice % (`.PermutationGroupRightCosets`): Right cosets in the permutation group
    end

    methods

        function self = RightCosets(group, subgroup)
            self.group = group;
            self.subgroup = subgroup;
            self.nice = group.niceGroup.rightCosetsOf(subgroup.niceGroup);
        end

        function t = canonicalRepresentative(self, g)
            t = self.group.niceMonomorphismPreimage(self.nice.canonicalRepresentative(self.group.niceMonomorphismImage(g)));
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
            T = cellfun(@(t) self.group.niceMonomorphismPreimage(t), self.nice.transversal, 'uniform', 0);
        end

    end

end
