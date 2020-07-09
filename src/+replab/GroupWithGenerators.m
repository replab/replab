classdef GroupWithGenerators < replab.Group
% Describes a group with generators
%
% It is used as a common base class between `.FiniteGroup` (which is compact), and
% `.FreeGroup` (which is not).
%
% The generators must never contain the identity.

    properties (SetAccess = protected)
        generators % row cell array of group elements: Group generators
    end

    methods

        function n = nGenerators(self)
        % Returns the number of group generators
        %
        % Returns:
        %   integer: Number of group generators
            n = length(self.generators);
        end

        function p = generator(self, i)
        % Returns the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: i-th group generator
            p = self.generators{i};
        end

        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: Inverse of the i-th group generator
            p = self.inverse(self.generators{i});
        end

        function b = isTrivial(self)
        % Tests whether this group is trivial
        %
        % Returns:
        %   logical: True if this group is trivial (i.e. has only one element)
            b = self.nGenerators == 0;
        end

    end

end
