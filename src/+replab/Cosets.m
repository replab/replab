classdef Cosets < replab.Obj

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        subgroup % (`.FiniteGroup`): Subgroup of `.group`
    end

end
