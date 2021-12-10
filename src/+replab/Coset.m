classdef Coset < replab.FiniteSet
% Base class for the left/normal/right coset classes
%
% This class is necessary because `.NormalCoset` inherits from both `.LeftCoset` or `.RightCoset`,
% and thus they share the `.group` and `.subgroup` properties.

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group containing this coset
        subgroup % (`.FiniteGroup`): Subgroup used to decompose `.group`
    end

end
