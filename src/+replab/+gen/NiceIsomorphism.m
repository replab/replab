classdef NiceIsomorphism < replab.FiniteIsomorphism
% An isomorphism from a subgroup of a finite group type to a finite group (Abstract)
%
% In addition, the isomorphism should preserve the type order.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.OrderPreservingFiniteIsomorphismLaws(self);
        end

    end

    methods

        function l = sourceContains(self, s)
        % Returns whether the source of this morphism contains the given type element
        %
        % Args:
        %   s (element of `.sourceType`): Element to check
        %
        % Returns:
        %   logical: True if source contains the given element
            error('Abstract');
        end

    end

end
