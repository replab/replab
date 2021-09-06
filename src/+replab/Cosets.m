classdef Cosets < replab.Obj

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        subgroup % (`.FiniteGroup`): Subgroup of `.group`
    end

    methods

        function self = Cosets(group, subgroup)
            assert(group.hasSameTypeAs(subgroup));
            self.group = group;
            self.subgroup = subgroup;
        end

    end

end
