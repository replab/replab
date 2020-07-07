classdef LeftCosets < replab.LeftCosets

    properties
        nice % (`.PermutationGroupLeftCosets`): Left cosets in the permutation group
    end

    methods

        function self = LeftCosets(group, subgroup)
            self.group = group;
            self.subgroup = subgroup;
            self.nice = group.niceGroup.leftCosetsOf(subgroup.niceGroup);
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
