classdef Cosets < replab.Obj

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        subgroup % (`.FiniteGroup`): Subgroup of `.group`
    end

    methods

        function s = nElements(self)
        % Returns the number of cosets
        %
        % Returns:
        %   integer: Number of cosets
            s = self.group.order / self.subgroup.order;
            assert(s <= 2^53 - 1);
            s = double(s);
        end

    end

end
