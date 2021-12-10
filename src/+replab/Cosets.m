classdef Cosets < replab.Obj
% Base class for the left/normal/right coset classes
%
% This class is necessary because `.NormalCosets` inherits from both `.LeftCosets` or `.RightCosets`,
% and thus they share the `.group` and `.subgroup` properties.

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
